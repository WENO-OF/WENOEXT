DOCKER_TAG_OF_v2406=wenotest:v2406

all: runTestsOFv2406


runTestsOFv2406:
	docker build \
		--rm \
		-t ${DOCKER_TAG_OF_v2406} \
		-f CI/Dockerfile.OFv2406 \
		.

removeStoppedContainer:
	docker container prune -f

clean: removeStoppedContainer
	docker rmi ${DOCKER_TAG_OF_v2406}
