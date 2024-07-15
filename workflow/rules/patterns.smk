# {input_dir}/{run}/{subrun}/{flowcell}/pod5
RGX_DIR_PATTERN = re.compile(
    (
        f'{config["input_dir"]}/'
        + r"(?P<run_dir>"
        + config["run_dir_pattern"]
        + r")/(?P<subrun_dir>[^/]+)/(?P<flowcell>[^/]+)/pod5$"
    )
)

OUTPUT_DIR = os.path.join(
    config["input_dir"],
    "{run_dir}",
    "{subrun_dir}",
    "{flowcell}",
    "pod5",
    config["run_output_dir"],
)

READS_FILE = os.path.join(
    OUTPUT_DIR,
    "{run_dir}.bam",
)
