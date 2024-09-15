def move_incomplete_files(
    incomplete_files: str | None, *, is_dry_run: bool
) -> dict[str, str]:
    # Exit if snakemake directory is locked
    try:
        if os.listdir(".snakemake/locks"):
            return {}
    except FileNotFoundError:
        pass

    # Check if need to rebasecall on failed Snakemake run.
    incomplete_file_paths = {}
    if not incomplete_files:
        return incomplete_file_paths

    with open(config["incomplete_files"], "rt") as fh:
        # Move incomplete files to /tmp. Snakemake deletes on new run otherwise.
        for l in fh.readlines():
            abs_path = l.strip()
            path_comp = abs_path.split("/")
            tmp_dir = os.path.join("/", "/".join(["tmp", *abs_path.split("/")[2:-1]]))
            new_abs_path = os.path.join(tmp_dir, path_comp[-1])
            os.makedirs(tmp_dir, exist_ok=True)
            incomplete_file_paths[abs_path] = new_abs_path
            # Move files.
            # If dry run or is not a bam file (basecalled file), continue.
            if is_dry_run or not abs_path.endswith(".bam"):
                continue
            try:
                shutil.move(abs_path, new_abs_path)
            except FileNotFoundError:
                pass

    print(
        "Incomplete files found:", list(incomplete_file_paths.keys()), file=sys.stderr
    )

    return incomplete_file_paths


def get_run_dir_wcs(rgx_dir_pattern: str) -> tuple[list[str], list[str], list[str]]:
    # Walk through input directory check that regex pattern matches
    run_dirs, subrun_dirs, flowcell_dirs = [], [], []
    for root, _, _ in os.walk(config["input_dir"]):
        mtch = re.search(rgx_dir_pattern, root)
        if not mtch:
            continue

        run_dir, subrun_dir, flowcell_dir = mtch.groups()

        abs_path_flowcell_dir = os.path.join(
            config["input_dir"], run_dir, subrun_dir, flowcell_dir
        )
        contains_summary_files = any(
            file.startswith("sequencing_summary")
            for file in os.listdir(abs_path_flowcell_dir)
        )
        # Skip directories that are still sequencing.
        if not contains_summary_files:
            print(
                f"No sequence summary file in {abs_path_flowcell_dir}. Skipping...",
                file=sys.stderr,
            )
            continue

        run_dirs.append(run_dir)
        subrun_dirs.append(subrun_dir)
        flowcell_dirs.append(flowcell_dir)

    return run_dirs, subrun_dirs, flowcell_dirs


def expand_path(wc, path: str):
    return expand(
        path, zip, run_dir=wc.run_dir, subrun_dir=wc.subrun_dir, flowcell=wc.flowcell
    )[0]
