import os
import shutil
import argparse
import logging
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def load_rename_map(tsv_file: str) -> dict:
    """
    Load two-column TSV into a { old_name: new_name } dict.
    Header row (any column order) is detected automatically.
    Without a header: first col = old, second col = new.
    """
    df = pd.read_csv(tsv_file, sep="\t", header=None, comment="#", dtype=str)

    first = df.iloc[0, 0].strip().lower()
    has_header = first in ("old_name", "old_names", "new_name", "new_names", "old", "new", "source")

    if has_header:
        header = [c.strip().lower().replace("_names", "_name") for c in df.iloc[0]]
        df = df.iloc[1:].reset_index(drop=True)
        df.columns = header
    else:
        df.columns = ["old_name", "new_name"]

    rename_map = dict(zip(df["old_name"].str.strip(), df["new_name"].str.strip()))

    if rename_map:
        sample_old, sample_new = next(iter(rename_map.items()))
        log.info(f"Loaded {len(rename_map)} entries. Sample: '{sample_old}' → '{sample_new}'")

    return rename_map


def apply_rename(name: str, rename_map: dict) -> str:
    """
    Replace old_name prefix in a file or folder name with new_name.

    Example:
        G472_FKDN260183859-1A_L1_1.fq.gz → AEO25_V2_251211_F26_FKDN260183859-1A_L1_1.fq.gz
        G472 (folder)                     → AEO25_V2_251211_F26
    """
    for old, new in rename_map.items():
        if old in name:
            parts = name.split(old, maxsplit=1)
            return new + parts[1]
    return None  # no match


def distribute(
    source: str,
    destination: str,
    rename_map: dict,
    dry_run: bool = False,
) -> None:
    """
    Walk source recursively. Only process top-level folders (G***) that are
    listed in rename_map — everything else is skipped.
    Files inside matched folders are renamed and copied to the destination.
    """
    skipped = []

    for entry in sorted(os.scandir(source), key=lambda e: e.name):
        if not entry.is_dir():
            continue

        folder = entry.name
        new_folder = apply_rename(folder, rename_map)

        if new_folder is None:
            skipped.append(folder)
            continue

        src_dir  = entry.path
        dest_dir = os.path.join(destination, new_folder)

        log.info(f"{'[DRY-RUN] ' if dry_run else ''}Dir: {folder} → {new_folder}")

        if not dry_run:
            os.makedirs(dest_dir, exist_ok=True)

        for filename in sorted(os.listdir(src_dir)):
            src_file  = os.path.join(src_dir, filename)
            new_fname = apply_rename(filename, rename_map) or filename
            dest_file = os.path.join(dest_dir, new_fname)

            if filename != new_fname:
                log.info(f"  Rename : {filename} → {new_fname}")
            else:
                log.info(f"  Copy   : {filename}")

            if dry_run:
                log.info(f"  [DRY-RUN] {src_file} → {dest_file}")
                continue

            if os.path.exists(dest_file):
                log.info(f"  Skipped (already exists): {dest_file}")
            else:
                shutil.copy2(src_file, dest_file)

    if skipped:
        log.info(f"Skipped {len(skipped)} folders not in TSV: {skipped}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="rename_distribute",
        description="Copy and rename a flat G*** folder structure via a TSV map, skipping unlisted folders.",
    )
    parser.add_argument("-s", "--source",      required=True, help="Source root directory (e.g. 01.RawData)")
    parser.add_argument("-d", "--destination", required=True, help="Destination root directory")
    parser.add_argument("-n", "--rename_tsv",  required=True, help="TSV with columns old_names / new_names (any order, header optional)")
    parser.add_argument("--dry_run", action="store_true",     help="Preview actions without copying")
    args = parser.parse_args()

    rename_map = load_rename_map(args.rename_tsv)

    distribute(
        source=args.source,
        destination=args.destination,
        rename_map=rename_map,
        dry_run=args.dry_run,
    )
