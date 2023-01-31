import os

import nox


@nox.session(python=['3.6', '3.9'], reuse_venv=False)
def run_tests(session):
    session.install('-r', 'requirements-dev.txt')
    session.run("pytest", '-n', '5')


@nox.session(python=False)
def serve_docs(session):
    session.run("sphinx-autobuild", "docs/source", "docs/build", "--port", "8777", "--open-browser")


@nox.session(python=False)
def build_docs(session):
    session.run("sphinx-build", "docs/source", "docs/build")


@nox.session(python=False)
def update_dev_kraken(session):
    """Run the Kraken build to update it with new pylinac changes"""
    key_info = os.environ['GOOGLE_CREDENTIALS']
    with open("service_key.json", "w") as key_file:
        key_file.write(key_info)
    session.run("gcloud", "auth", "activate-service-account", "--key-file service_key.json")
    session.run("gcloud", "config", "set", "project", 'radmachine')
    session.run("gcloud", "builds", "triggers", "run", "kraken-build", "--branch=master")
