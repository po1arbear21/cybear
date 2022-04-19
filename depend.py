import os, argparse, collections, sys, re

# regex programs
re_include   = re.compile(r"m4_include\(\s*(.*)\)")
re_program   = re.compile(r"^program *(.*?)(?: *!|$)")
re_module    = re.compile(r"^module *(.*?)(?: *!|$)")
re_submodule = re.compile(r"^submodule *\((.*)\) *(.*?)(?: *!|$)")
re_use       = re.compile(r"^use(?:.*::)? *(.*?)(?: *,|$| *!)")
re_library   = re.compile(r"^\!\!\!library *(.*)")

class ParseException(Exception):
  def __init__(self, file, line_number, line, msg):
    self.file        = file
    self.line_number = line_number
    self.line        = line
    self.msg         = msg

class src_item:
  def __init__(self, name, file_index):
    self.name       = name
    self.file_index = file_index

class m4_file:
  def __init__(self, folder, filename):
    self.folder   = folder
    self.filename = filename
    self.includes = []

    file, extension = os.path.splitext(filename)
    if (extension != ".f90"): raise Exception("file extension must be *.f90")

    full_filename = folder + "/" + filename

    # open file
    with open(full_filename) as f:
      # start at the top
      line_number = 0

      # read and parse lines to get includes
      while True:
        line = f.readline()
        if (line == ""): break
        line_number = line_number + 1

        # remove leading and trailing spaces
        line = line.strip(" \r\n")

        # search for m4_include
        result = re_include.match(line)
        if (result is None): continue

        # add included file to includes list
        inc = result.group(1)
        if (inc == ""): raise ParseException(full_filename, line_number, line, "m4_include filename empty")
        if (inc not in self.includes): self.includes.append(inc)

class f90_file:
  def __init__(self, folder, filename):
    self.folder       = folder
    self.filename     = filename
    self.program_name = ""
    self.library_name = ""
    self.modules      = []
    self.submodules   = []
    self.use_modules  = []
    self.dep_anchor   = []
    self.parents      = []

    # object and anchor name
    file, extension = os.path.splitext(filename)
    if (extension != ".f90"): raise Exception("filename must be *.f90")
    self.objname = file + ".o"
    self.ancname = file + ".anc"

    full_filename = folder + "/" + filename

    # open file
    with open(full_filename) as f:
      # start at the top
      line_number = 0

      # read and parse lines
      while True:
        line = f.readline()
        if (line == ""): break
        line_number = line_number + 1

        # remove leading and trailing spaces
        line = line.strip(" \r\n")

        # check for library statement
        result = re_library.match(line)
        if (result is not None):
          if (result.group(1) == ""): raise ParseException(full_filename, line_number, line, "library name missing")
          if (self.library_name != ""): raise ParseException(full_filename, line_number, line, "multiple library statements in one file")
          if (self.program_name != ""): raise ParseException(full_filename, line_number, line, "program and library statement in one file")
          self.library_name = result.group(1).lower()

          continue

        # check for program statement
        result = re_program.match(line)
        if (result is not None):
          if (result.group(1) == ""): raise ParseException(full_filename, line_number, line, "program name missing")
          if (self.program_name != ""): raise ParseException(full_filename, line_number, line, "multiple program statements in one file")
          if (self.library_name != ""): raise ParseException(full_filename, line_number, line, "program and library statement in one file")
          self.program_name = result.group(1).lower()

          continue

        # check for module statement
        result = re_module.match(line)
        if (result is not None):
          if (result.group(1) == ""): raise ParseException(full_filename, line_number, line, "module name missing")
          module_name = result.group(1).lower()

          # make sure this is not a "module subroutine ..." etc. statement
          if (module_name.find(" ") == -1):
            if (module_name not in self.modules): self.modules.append(module_name)

          continue

        # check for submodule statement
        result = re_submodule.match(line)
        if (result is not None):
          if (result.group(1) == ""): raise ParseException(full_filename, line_number, line, "submodule parent missing")
          parent = result.group(1).lower()
          if (result.group(2) == ""): raise ParseException(full_filename, line_number, line, "submodule name missing")
          submodule_name = result.group(2).lower()
          submodule_name = parent + "@" + submodule_name
          if (parent not in self.use_modules): self.use_modules.append(parent)
          if (parent not in self.parents): self.parents.append(parent)
          if (submodule_name not in self.submodules): self.submodules.append(submodule_name)

          continue

        # check for use statement
        result = re_use.match(line)
        if (result is not None):
          if (result.group(1) == ""): raise ParseException(full_filename, line_number, line, "use module name missing")
          module_name = result.group(1).lower()
          if (module_name not in self.use_modules): self.use_modules.append(module_name)

          continue

    # init module dependency index list
    self.dep_modules = [-1] * len(self.use_modules)

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("filenames", nargs = "+")
args = parser.parse_args()

files = []
for f in args.filenames:
  folder, filename = os.path.split(f)

  # m4
  file = m4_file(folder, filename)
  print("$(BUILD_DIR)/" + filename + ": " + folder + "/" + filename, end = "")
  for i in range(len(file.includes)):
    if (os.path.isfile(folder + "/" + file.includes[i])):
      print (" " + folder + "/" + file.includes[i], end = "")
  print("")

  # collect f90 files
  files.append(f90_file(folder, filename))

# collect all programs, libraries, modules and submodules
programs   = []
libraries  = []
modules    = []
submodules = []
for i in range(len(files)):
  if (files[i].program_name != ""):
    for j in range(len(programs)):
      if (files[i].program_name == programs[j].name): raise Exception("program name " + programs[j].name + " found more than once")
    programs.append(src_item(files[i].program_name, i))

  if (files[i].library_name != ""):
    for j in range(len(libraries)):
      if (files[i].library_name == libraries[j].name): raise Exception("library name " + libraries[j].name + " found more than once")
    libraries.append(src_item(files[i].library_name, i))

  for j in range(len(files[i].modules)):
    for k in range(len(modules)):
      if (modules[k].name == files[i].modules[j]): raise Exception("module name " + modules[k].name + " found more than once")
    modules.append(src_item(files[i].modules[j], i))

  for j in range(len(files[i].submodules)):
    for k in range(len(submodules)):
      if (submodules[k].name == files[i].submodules[j]): raise Exception("submodule name " + submodules[k].name + " found more than once")
    submodules.append(src_item(files[i].submodules[j], i))

# get dependency indices
for i in range(len(files)):
  for j in range(len(files[i].use_modules)):
    for k in range(len(modules)):
      if (modules[k].name == files[i].use_modules[j]):
        files[i].dep_modules[j] = k

        # add anchor file of this module (if it does not exist already)
        if (modules[k].file_index not in files[i].dep_anchor):
          files[i].dep_anchor.append(modules[k].file_index)

        break

# get object list for each program
pobjects = []
for i in range(len(programs)):
  pobjects.append([programs[i].file_index])

  # init deque of file candidates
  file_deque = collections.deque()
  for d in files[programs[i].file_index].dep_anchor:
    file_deque.append(d)

  # add files and their dependencies
  while (len(file_deque) > 0):
    j = file_deque.popleft()

    # add file if not already in list
    if (j not in pobjects[i]): pobjects[i].append(j)

    # add dependencies
    for d in files[j].dep_anchor:
      file_deque.append(d)

  # add submodules
  for s in submodules:
    for j in range(len(files[s.file_index].parents)):
      for k in pobjects[i]:
        if (files[s.file_index].parents[j] in files[k].modules):
          if (s.file_index not in pobjects[i]): pobjects[i].append(s.file_index)

# get object list for each library
lobjects = []
for i in range(len(libraries)):
  lobjects.append([libraries[i].file_index])

  # init deque of file candidates
  file_deque = collections.deque()
  for d in files[libraries[i].file_index].dep_anchor:
    file_deque.append(d)

  # add files and their dependencies
  while (len(file_deque) > 0):
    j = file_deque.popleft()

    # add file if not already in list
    if (j not in lobjects[i]): lobjects[i].append(j)

    # add dependencies
    for d in files[j].dep_anchor:
      file_deque.append(d)

  # add submodules
  for s in submodules:
    for j in range(len(files[s.file_index].parents)):
      for k in lobjects[i]:
        if (files[s.file_index].parents[j] in files[k].modules):
          if (s.file_index not in lobjects[i]): lobjects[i].append(s.file_index)

# output anchor files
print("")
for i in range(len(files)):
  "$(BUILD_DIR)/"
  print("$(BUILD_DIR)/" + files[i].ancname + ": $(BUILD_DIR)/" + files[i].filename, end = "")

  # depend on anchor files
  for j in range(len(files[i].dep_anchor)):
    print(" $(BUILD_DIR)/" + files[files[i].dep_anchor[j]].ancname, end = "")

  # next line
  print("")

# output object files
print("")
for i in range(len(files)):
  print("$(BUILD_DIR)/" + files[i].objname + ": $(BUILD_DIR)/" + files[i].filename, end = "")

  # depend on anchor files
  for j in range(len(files[i].dep_anchor)):
    print(" $(BUILD_DIR)/" + files[files[i].dep_anchor[j]].ancname, end = "")

  # next line
  print("")

# output anchor list for each target
print("")
for i in range(len(programs)):
  print("ANCHORS_" + programs[i].name + " =", end = "")
  for j in pobjects[i]:
    print(" $(BUILD_DIR)/" + files[j].ancname, end = "")
  print("")
for i in range(len(libraries)):
  print("ANCHORS_" + libraries[i].name + " =", end = "")
  for j in lobjects[i]:
    print(" $(BUILD_DIR)/" + files[j].ancname, end = "")
  print("")

# output object list for each target
print("")
for i in range(len(programs)):
  print("OBJECTS_" + programs[i].name + " =", end = "")
  for j in pobjects[i]:
    print(" $(BUILD_DIR)/" + files[j].objname, end = "")
  print("")
for i in range(len(libraries)):
  print("OBJECTS_" + libraries[i].name + " =", end = "")
  for j in lobjects[i]:
    print(" $(BUILD_DIR)/" + files[j].objname, end = "")
  print("")

# output targets variable (programs and libraries)
print("")
print("TARGETS =", end = "")
for i in range(len(programs)):
  print(" $(BUILD_DIR)/" + programs[i].name, end = "")
for i in range(len(libraries)):
  print(" $(BUILD_DIR)/" + libraries[i].name, end = "")

# output targets
print("\n")
for i in range(len(programs)):
  print(".PHONY: " + programs[i].name)
  print("$(BUILD_DIR)/" + programs[i].name + ": $(ANCHORS_" + programs[i].name + ") $(OBJECTS_" + programs[i].name + ") $(C_OBJECTS)")
  print('\t@printf "%b" "$(FC_P) $(FFLAGS) $(FINCLUDE) -o $(OUT_COLOR)' + "$(BUILD_DIR)/" + programs[i].name + '$(OFF_COLOR) $(INP_COLOR)$(OBJECTS_' + programs[i].name + ') $(C_OBJECTS) $(LIBS)$(OFF_COLOR)\\n\\n"')
  print("\t@$(FC) $(FFLAGS) $(FINCLUDE) -o $(BUILD_DIR)/" + programs[i].name + " $(OBJECTS_" + programs[i].name + ") $(C_OBJECTS) $(LIBS)")
  print("")
for i in range(len(libraries)):
  print(".PHONY: " + libraries[i].name)
  print("$(BUILD_DIR)/" + libraries[i].name + ": $(ANCHORS_" + libraries[i].name + ") $(OBJECTS_" + libraries[i].name + ") $(C_OBJECTS)")
  print('\t@printf "%b" "$(AR_P) rcs $(OUT_COLOR)' + "$(BUILD_DIR)/" + libraries[i].name + '.a$(OFF_COLOR) $(INP_COLOR)$(OBJECTS_' + libraries[i].name + ') $(C_OBJECTS)$(OFF_COLOR)\\n\\n"')
  print("\t@ar rcs $(BUILD_DIR)/" + libraries[i].name + ".a $(OBJECTS_" + libraries[i].name + ") $(C_OBJECTS)")
  print("")
