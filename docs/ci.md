# CI

## Manage the fargo docker image

### Login to registry.git.rwth-aachen.de

Go to [the project settings](https://git.rwth-aachen.de/ithe/fargo/-/settings/access_tokens)
for fargo and create an access token with `{read,write}_registry` permissions.
Login with docker locally:

```bash
echo <token> | docker login registry.git.rwth-aachen.de -u <gitlab-user> --password-stdin
```

### Create an image

```bash
docker buildx build -t registry.git.rwth-aachen/ithe/fargo
```

### Upload the image

```bash
docker push registry.git.rwth-aachen/ithe/fargo
```

## Using the image with GitLab CI/CD

## Getting registry access

By default, other projects don't have access to a project's container registry.
To grant access for fortran-basic, [the job token permission](https://git.rwth-aachen.de/ithe/fargo/-/settings/ci_cd#js-token-access)
page for fargo and add fortran-basic to the job token allowlist.

## Example usage

Each CI run generates a short-lived job token which will be used to log in to
private GitLab instances like git.rwth-aachen.de and grant access to several
features including the docker registries. Since the project was granted access
in the previous step, the image can simply be used inside the `.gitlab-ci.yml` file.

```yaml
...
test:
  image: registry.git.rwth-aachen.de/ithe/fargo

...
```
