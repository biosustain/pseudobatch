.PHONY: qa

## Apply code quality assurance tools.
qa:
	isort --recursive pseudobatch/ tests/
	black pseudobatch/ tests/

