DOCKER_TAG_OF_5=wenotest:5.x
DOCKER_TAG_OF_7=wenotest:7
DOCKER_TAG_OF_8=wenotest:8
DOCKER_TAG_OF_v1912=openfoam:v1912

all: runTestsOF5 runTestsOF7 runTestsOF8


runTestsOF5:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_5} \
		-f CI/Dockerfile.OF5 \
		.

runTestsOF7:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_7} \
		-f CI/Dockerfile.OF7 \
		.

runTestsOF8:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_8} \
		-f CI/Dockerfile.OF8 \
		.

removeStoppedContainer:
	docker container prune -f

clean: removeStoppedContainer
	docker rmi ${DOCKER_TAG_OF_5} ${DOCKER_TAG_OF_7} $ ${DOCKER_TAG_OF_8}
