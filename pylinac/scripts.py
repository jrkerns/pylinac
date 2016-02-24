"""Scripts to run pylinac from the command line."""
import os.path as osp

from pylinac.watcher import start_watching
try:
    import click
except ImportError as e:
    e.msg = "To run pylinac scripts you must have `click` installed. Run `pip install click`, then try again."
    raise e


@click.group()
def cli():
    """Perform pylinac operations from the command line."""
    pass


@click.command()
@click.argument('directory', type=click.Path(exists=True))
@click.option('--config', type=click.Path(exists=True))
def watch(directory, config=None):
    """Start watching a directory and analyze any applicable files"""
    start_watching(directory, config)


cli.add_command(watch)

if __name__ == '__main__':
    the_dir = osp.dirname(__file__)
    watch((the_dir,))
