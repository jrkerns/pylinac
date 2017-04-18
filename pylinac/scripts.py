"""Scripts to run pylinac from the command line. Built on examples given in Click documentation."""
from pylinac import watcher
from pylinac import log_analyzer
from pylinac import gui as pgui
import click


@click.group()
def cli():
    """Perform pylinac operations from the command line."""
    pass


@cli.command()
@click.option('--directory', type=click.Path(exists=True))
@click.option('--config', type=click.Path(exists=True))
def watch(directory=None, config=None):
    """Start watching a directory and analyze any applicable files"""
    watcher.watch(directory, config)


@cli.command()
@click.option('--directory', type=click.Path(exists=True))
@click.option('--config', type=click.Path(exists=True))
@click.option('--move-new-logs', type=click.BOOL)
def process(directory=None, config=None):
    """Process any applicable files in the directory once through."""
    watcher.process(directory, config)


@cli.command()
@click.argument('directory', type=click.Path(exists=True))
@click.option('--destination', type=click.Path(exists=True), help='Destination folder to place anonymized logs.')
@click.option('--in-place', is_flag=True, help='Whether to modify logs in-place or create copies.')
def anonymize(directory, destination=None, in_place=None):
    """Anonymize machine logs in a given file directory."""
    log_analyzer.anonymize(directory, inplace=in_place, destination=destination)


@cli.command()
def gui():
    """Start the Pylinac GUI."""
    pgui()
