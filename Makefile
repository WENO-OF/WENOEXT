DOCKER_TAG_OF_5=wenotest:5.x
DOCKER_TAG_OF_7=wenotest:7
DOCKER_TAG_OF_8=wenotest:8
DOCKER_TAG_OF_v1912=wenotest:v1912
DOCKER_TAG_OF_v2006=openfoam:v2006weno
DOCKER_TAG_OF_v2012=openfoam:v2012weno
DOCKER_TAG_OF_v2506=wenotest:v2506
DOCKER_TAG_OF_v2512=wenotest:v2512

all: runTestsOF5 runTestsOF7 runTestsOF8 runTestsOFv1912 runTestsOFv2006 runTestsOFv2012 runTestsOFv2506 runTestsOFv2512


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

runTestsOFv1912:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_v1912} \
		-f CI/Dockerfile.OFv1912 \
		.

runTestsOFv2006:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_v2006} \
		-f CI/Dockerfile.OFv2006 \
		.
runTestsOFv2012:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_v2012} \
		-f CI/Dockerfile.OFv2012 \
		.

runTestsOFv2506:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_v2506} \
		-f CI/Dockerfile.OFv2506 \
		.

runTestsOFv2512:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_v2512} \
		-f CI/Dockerfile.OFv2512 \
		.

removeStoppedContainer:
	docker container prune -f

clean: removeStoppedContainer
	docker rmi ${DOCKER_TAG_OF_5} ${DOCKER_TAG_OF_7} ${DOCKER_TAG_OF_8} ${DOCKER_TAG_OF_v1912} ${DOCKER_TAG_OF_v2006} ${DOCKER_TAG_OF_v2012} ${DOCKER_TAG_OF_v2506} ${DOCKER_TAG_OF_v2512}
