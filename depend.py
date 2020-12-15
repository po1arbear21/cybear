import os, argparse, queue, sys

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

class src_file:
  def __init__(self, folder, filename):
    self.folder       = folder
    self.filename     = filename
    self.program_name = ""
    self.library_name = ""
    self.modules      = []
    self.submodules   = []
    self.use_modules  = []
    self.includes     = []
    self.dep_anchor   = []
    self.parents      = []

    # object and anchor name
    file, extension = os.path.splitext(filename)
    if (extension != ".f90"): raise Exception("filename must be *.f90")
    self.objname = file + ".o"
    self.ancname = file + ".anc"

    full_filename = folder + filename

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

        # check for library statement
        if ((len(line) >= 11) and (line[0:11].lower() == "!!!library ")):
          if (len(line) < 12): raise ParseException(full_filename, line_number, line, "library name missing")
          if (self.library_name != ""): raise ParseException(full_filename, line_number, line, "multiple library statements in one file")
          if (self.program_name != ""): raise ParseException(full_filename, line_number, line, "program and library statement in one file")
          self.library_name = line[11:].lower().strip()

          continue

        # remove comments
        quote = ""
        for i in range(len(line)):
          if (quote == ""): # not inside quoted string
            # check for beginning of quoted string or comment
            if ((line[i] == "'") or (line[i] == '"')):
              # save quote character
              quote = line[i]
            elif (line[i] == "!"):
              # trim line up to comment character
              line = line[0:i-1]
              break
          elif (line[i] == quote):
            # end of quoted string
            quote = ""

        # check for program statement
        if ((len(line) >= 8) and (line[0:8].lower() == "program ")):
          if (len(line) < 9): raise ParseException(full_filename, line_number, line, "program name missing")
          if (self.program_name != ""): raise ParseException(full_filename, line_number, line, "multiple program statements in one file")
          if (self.library_name != ""): raise ParseException(full_filename, line_number, line, "program and library statement in one file")
          self.program_name = line[8:].lower().strip()

          continue

        # check for module statement
        if ((len(line) >= 7) and (line[0:7].lower() == "module ")):
          if (len(line) < 8): raise ParseException(full_filename, line_number, line, "module name missing")
          module_name = line[7:].strip().lower()

          # make sure this is not a "module procedure ..." etc. statement
          if (module_name.find(" ") == -1):
            if (module_name not in self.modules): self.modules.append(module_name)

          continue

        # check for submodule statement
        if ((len(line) >= 10) and (line[0:9].lower() == "submodule") and ((line[9] == " ") or line[9] == "(")):
          # get rid of "submodule"
          submodule_name = line[9:].strip()

          # get parent name
          if (submodule_name[0] != "("): raise ParseException(full_filename, line_number, line, "submodule parent missing")
          i = submodule_name.find(")")
          if (i <= 1): raise ParseException(full_filename, line_number, line, "closing parenthesis missing")
          parent = submodule_name[1:i].strip()
          if (parent == ""): raise ParseException(full_filename, line_number, line, "submodule parent empty")
          parent = parent.lower()

          # get submodule name
          if (len(submodule_name) <= i + 1): raise ParseException(full_filename, line_number, line, "submodule name missing")
          submodule_name = line[i+1:].strip()
          if (submodule_name == ""): raise ParseException(full_filename, line_number, line, "submodule name missing")
          submodule_name = parent + "@" + submodule_name.lower()

          # use parent module
          if (parent not in self.use_modules): self.use_modules.append(parent)
          if (parent not in self.parents): self.parents.append(parent)

          # add submodule
          if (submodule_name not in self.submodules): self.submodules.append(submodule_name)

          continue

        # check for use statement
        if ((len(line) >= 4) and (line[0:4] == "use ")):
          if (len(line) < 5): raise ParseException(full_filename, line_number, line, "use module name missing")

          # get rid of "use "
          module_name = line[4:].strip()
          if (module_name == ""): raise ParseException(full_filename, line_number, line, "use module name missing")

          # search and remove ::
          i = module_name.find("::")
          if (len(module_name) <= i + 2): raise ParseException(full_filename, line_number, line, "use module name missing")
          if (i != -1): module_name = module_name[i+2:].strip()

          # remove stuff after ,
          i = module_name.find(",")
          if (i != -1):
            if (i < 1): raise ParseException(full_filename, line_number, line, "use module name missing")
            module_name = module_name[0:i].strip()

          if (module_name == ""): raise ParseException(full_filename, line_number, line, "use module name missing")
          module_name = module_name.lower()

          if (module_name not in self.use_modules): self.use_modules.append(module_name)

          continue

        # check for includes
        if ((len(line) >= 7) and (line[0:7] == "include")):
          inc = line[7:].strip()
        elif ((len(line) >= 8) and (line[0:8] == "#include")):
          inc = line[8:].strip()
        else:
          inc = ""
        if (inc != ""):
          if (len(inc) < 3): raise ParseException(full_filename, line_number, line, "include filename missing")
          quote = inc[0]
          if ((quote != "'") and (quote != '"')): raise ParseException(full_filename, line_number, line, "include filename quote missing")
          if (inc[-1] != quote): raise ParseException(full_filename, line_number, line, "include filename closing quote missing")
          inc = inc[1:-1].strip()
          if (inc == ""): raise ParseException(full_filename, line_number, line, "include filename empty")
          if (inc not in self.includes): self.includes.append(inc)

          continue

    # init module dependency index list
    self.dep_modules = [-1] * len(self.use_modules)

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--build_dir", required = True)
parser.add_argument("-I", "--include_dirs", nargs = "+")
parser.add_argument("filenames", nargs = "+")
args = parser.parse_args()

# initialize source files
files = []
for f in args.filenames:
  folder, filename = os.path.split(f)
  files.append(src_file(folder + "/", filename))

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

  # search include files in include folders
  for j in range(len(files[i].includes)):
    file_exists = False

    # search file folder first
    if (os.path.isfile(files[i].folder + files[i].includes[j])):
      files[i].includes[j] = files[i].folder + files[i].includes[j]
      file_exists = True
    else:
      if (args.include_dirs != None):
        # search include folders in order
        for k in range(len(args.include_dirs)):
          if (os.path.isfile(args.include_dirs[k] + files[i].includes[j])):
            files[i].includes[j] = args.include_dirs[k] + files[i].includes[j]
            file_exists = True
            break

    # do not include file if not found in folders (e.g. includes of a library interface)
    if (not file_exists): files[i].includes[j] = ""

# get object list for each program
pobjects = []
for i in range(len(programs)):
  pobjects.append([programs[i].file_index])

  # init queue of file candidates
  file_queue = queue.Queue()
  for d in files[programs[i].file_index].dep_anchor:
    file_queue.put(d)

  # add files and their dependencies
  while (not file_queue.empty()):
    j = file_queue.get()

    # add file if not already in list
    if (j not in pobjects[i]): pobjects[i].append(j)

    # add dependencies
    for d in files[j].dep_anchor:
      file_queue.put(d)

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

  # init queue of file candidates
  file_queue = queue.Queue()
  for d in files[libraries[i].file_index].dep_anchor:
    file_queue.put(d)

  # add files and their dependencies
  while (not file_queue.empty()):
    j = file_queue.get()

    # add file if not already in list
    if (j not in lobjects[i]): lobjects[i].append(j)

    # add dependencies
    for d in files[j].dep_anchor:
      file_queue.put(d)

  # add submodules
  for s in submodules:
    for j in range(len(files[s.file_index].parents)):
      for k in lobjects[i]:
        if (files[s.file_index].parents[j] in files[k].modules):
          if (s.file_index not in lobjects[i]): lobjects[i].append(s.file_index)

# output targets variable (programs and libraries)
print("TARGETS =", end = "")
for i in range(len(programs)):
  print(" " + args.build_dir + programs[i].name, end = "")
for i in range(len(libraries)):
  print(" " + args.build_dir + libraries[i].name, end = "")
print("\n")

# output object list for each target
for i in range(len(programs)):
  print("OBJECTS_" + programs[i].name + " =", end = "")
  for j in pobjects[i]:
    print(" " + args.build_dir + files[j].objname, end = "")
  print("")
for i in range(len(libraries)):
  print("OBJECTS_" + libraries[i].name + " =", end = "")
  for j in lobjects[i]:
    print(" " + args.build_dir + files[j].objname, end = "")
  print("")
print("")

# output anchor list for each target
for i in range(len(programs)):
  print("ANCHORS_" + programs[i].name + " =", end = "")
  for j in pobjects[i]:
    print(" " + args.build_dir + files[j].ancname, end = "")
  print("")
for i in range(len(libraries)):
  print("ANCHORS_" + libraries[i].name + " =", end = "")
  for j in lobjects[i]:
    print(" " + args.build_dir + files[j].ancname, end = "")
  print("")
print("")

# output targets
for i in range(len(programs)):
  print(".PHONY: " + programs[i].name)
  print(args.build_dir + programs[i].name + ": $(ANCHORS_" + programs[i].name + ") $(OBJECTS_" + programs[i].name + ") $(OBJECTS_C)")
  print('\t@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FFLAGS) $(FINT64) $(FREAL64) $(FINCLUDE) -o $(OU_COL)' + args.build_dir + programs[i].name + '$(NO_COL) $(IN_COL)$(OBJECTS_' + programs[i].name + ') $(OBJECTS_C) $(LIBS)$(NO_COL)\\n\\n"')
  print("\t@$(FC) $(FFLAGS) $(FINT64) $(FREAL64) $(FINCLUDE) -o " + args.build_dir + programs[i].name + " $(OBJECTS_" + programs[i].name + ") $(OBJECTS_C) $(LIBS)")
  print("")
for i in range(len(libraries)):
  print(".PHONY: " + libraries[i].name)
  print(args.build_dir + libraries[i].name + ": $(ANCHORS_" + libraries[i].name + ") $(OBJECTS_" + libraries[i].name + ") $(OBJECTS_C)")
  print('\t@printf "%b" "$(FC_COL)ar$(NO_COL) rcs $(OU_COL)' + args.build_dir + libraries[i].name + '.a$(NO_COL) $(IN_COL)$(OBJECTS_' + libraries[i].name + ') $(OBJECTS_C)$(NO_COL)\\n\\n"')
  print("\t@ar rcs " + args.build_dir + libraries[i].name + ".a $(OBJECTS_" + libraries[i].name + ") $(OBJECTS_C)")
  print("")

# output anchor files
for i in range(len(files)):
  print(args.build_dir + files[i].ancname + ": " + files[i].folder + files[i].filename, end = "")

  # depend on include files
  for j in range(len(files[i].includes)):
    if (files[i].includes[j] != ""):
      print(" " + files[i].includes[j], end = "")

  # depend on anchor files
  for j in range(len(files[i].dep_anchor)):
    print(" " + args.build_dir + files[files[i].dep_anchor[j]].ancname, end = "")

  # next line
  print("")

# output object files
for i in range(len(files)):
  print(args.build_dir + files[i].objname + ": " + files[i].folder + files[i].filename, end = "")

  # depend on include files
  for j in range(len(files[i].includes)):
    if (files[i].includes[j] != ""):
      print(" " + files[i].includes[j], end = "")

  # depend on anchor files
  for j in range(len(files[i].dep_anchor)):
    print(" " + args.build_dir + files[files[i].dep_anchor[j]].ancname, end = "")

  # next line
  print("")
