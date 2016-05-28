"""Scripts to run pylinac from the command line. Built on examples given in Click documentation."""
import argparse
import os.path as osp
import sys

from pylinac.watcher import start_watching
from pylinac import log_analyzer
try:
    import click
    from click.testing import CliRunner
except ImportError as e:
    e.msg = "To run pylinac scripts you must have `click` installed. Run `pip install click`, then try again."
    raise e


@click.group()
def cli():
    """Perform pylinac operations from the command line."""
    pass


@cli.command()
@click.argument('directory', type=click.Path(exists=True))
@click.option('--config', type=click.Path(exists=True))
def watch(directory, config=None):
    """Start watching a directory and analyze any applicable files"""
    start_watching(directory, config)


@cli.command()
@click.argument('directory', type=click.Path(exists=True))
@click.option('--destination', type=click.Path(exists=True), help='Destination folder to place anonymized logs.')
@click.option('--in-place', is_flag=True, help='Whether to modify logs in-place or create copies.')
def anonymize(directory, destination=None, in_place=None):
    """Anonymize machine logs in a given file directory."""
    log_analyzer.anonymize(directory, inplace=in_place, destination=destination)


if __name__ == '__main__':
    try:
        watch((sys.argv[1],))
    except IndexError:
        watch((osp.dirname(__file__),))
