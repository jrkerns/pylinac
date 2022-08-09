import nox


@nox.session(python=['3.6', '3.9'], reuse_venv=False)
def run_tests(session):
    session.install('-r', 'requirements-dev.txt')
    session.run("pytest", '-n', '5')


@nox.session(python=False)
def serve_docs(session):
    session.run("sphinx-autobuild", "docs/source", "docs/build", "--port", "8777")


@nox.session(python=False)
def build_docs(session):
    session.run("sphinx-build", "docs/source", "docs/build")
