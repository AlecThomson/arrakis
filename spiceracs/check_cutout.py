from tqdm import tqdm
import argparse
from glob import glob

def main(total, checkfile):
    with tqdm(total=total) as pbar:
        while True:
            current = len(glob(f'{checkfile}/*/*.fits'))
            pbar.update(current - pbar.n)
            if current == total:
                break

def cli():
        # Help string to be shown using the -h option
    logostr = """
     mmm   mmm   mmm   mmm   mmm
     )-(   )-(   )-(   )-(   )-(
    ( S ) ( P ) ( I ) ( C ) ( E )
    |   | |   | |   | |   | |   |
    |___| |___| |___| |___| |___|
     mmm     mmm     mmm     mmm
     )-(     )-(     )-(     )-(
    ( R )   ( A )   ( C )   ( S )
    |   |   |   |   |   |   |   |
    |___|   |___|   |___|   |___|

    """

    descStr = f"""
    {logostr}

    Monitor cutout progress

    """
    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'total',
        metavar='total',
        type=float,
        help='Total number of cutouts.')

    parser.add_argument(
        'checkfile',
        metavar='checkfile',
        type=str,
        help='Directory to monitor.')

    args = parser.parse_args()

    main(args.total, args.checkfile)

if __name__ == "__main__":
    cli()