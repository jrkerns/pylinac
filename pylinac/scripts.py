"""Scripts to run pylinac from the command line."""
import os.path as osp

from pylinac.watcher import start_watching
try:
    import click
except ImportError as e:
    raise e("To run scripts you must have click installed. Run `pip install click`, then try again.")

ANALYSIS_OPTIONS = ['cbct', 'vmat', 'starshot', 'leeds', 'pipspro-qc3', 'picketfence', 'log', 'winston-lutz']


@click.group()
def cli():
    """Perform pylinac operations from the command line."""
    pass


@click.command()
@click.argument('directory', type=click.Path(exists=True))
def watch(directory):
    """Start watching a directory and analyze any applicable files"""
    start_watching(directory)


@click.command()
@click.argument('target', type=click.Path(exists=True))
@click.option('--analysis-type', type=click.Choice(ANALYSIS_OPTIONS))
@click.option('--zip', is_flag=True)
def analyze(target, zip):
    pass


cli.add_command(watch)
cli.add_command(analyze)

if __name__ == '__main__':
    the_dir = osp.dirname(__file__)
    # watch((the_dir,))
    watch('_')
