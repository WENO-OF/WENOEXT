# Dockerfiles for Contineous Integration

Dockerfiles to generate the dockerimage to execute the unit tests with,

 * OpenFOAM v5.x
 * OpenFOAM v7
 * OpenFOAM v8
 * OpenFOAM v1912
 * OpenFOAM v2012


## Error Handling

If you cannot pull the required docker images for testing you can also
build the docker image yourself by using the provided template Dockerfile.

For generating ESI OpenFOAM docker images adapt the OpenFOAM version in the 
Dockerfile.createOFContainer file. For non ESI images the git clone paths 
need also to be updated.

The created docker images are also suitable for GitLab CI which is why a gitlab user
and goup are created. 
