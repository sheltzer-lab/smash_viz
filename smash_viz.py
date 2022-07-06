import argparse
from itertools import product
import pyparsing as pp
import re

# In order to highlight by arm we need to know the relative length of each arm.
_arms = {
    "1": {"p": 0.5015, "q": 0.4985},
    "2": {"p": 0.3836, "q": 0.6164},
    "3": {"p": 0.4595, "q": 0.5405},
    "4": {"p": 0.2637, "q": 0.7363},
    "5": {"p": 0.2675, "q": 0.7325},
    "6": {"p": 0.3565, "q": 0.6435},
    "7": {"p": 0.3764, "q": 0.6236},
    "8": {"p": 0.3116, "q": 0.6884},
    "9": {"p": 0.3470, "q": 0.6530},
    "10": {"p": 0.2966, "q": 0.7034},
    "11": {"p": 0.3978, "q": 0.6022},
    "12": {"p": 0.2675, "q": 0.7325},
    "13": {"p": 0.1554, "q": 0.8446},
    "14": {"p": 0.1640, "q": 0.8360},
    "15": {"p": 0.1853, "q": 0.8147},
    "16": {"p": 0.4051, "q": 0.5949},
    "17": {"p": 0.2956, "q": 0.7044},
    "18": {"p": 0.2203, "q": 0.7797},
    "19": {"p": 0.4482, "q": 0.5518},
    "20": {"p": 0.4363, "q": 0.5637},
    "21": {"p": 0.2743, "q": 0.7257},
    "22": {"p": 0.2865, "q": 0.7135},
    "X": {"p": 0.3903, "q": 0.6097},
    "Y": {"p": 0.2105, "q": 0.7895},
}

_rgb = {
    "black": (0, 0, 0),
    "white": (255, 255, 255),
    "red": (255, 0, 0),
    "lime": (0, 255, 0),
    "blue": (0, 0, 255),
    "yellow": (255, 255, 0),
    "cyan": (0, 255, 255),
    "magenta": (255, 0, 255),
    "silver": (192, 192, 192),
    "gray": (128, 128, 128),
    "maroon": (128, 0, 0),
    "olive": (128, 128, 0),
    "green": (0, 128, 0),
    "purple": (128, 0, 128),
    "teal": (0, 128, 128),
    "navy": (0, 0, 128),
    "lightblue": (173, 216, 230),
}

# Parser for chromosome lines
# Example:
# 64.909 47.242 m glx sp
# 92.805 564.9 m (1) jc s
# 120.7 47.242 m glx sp
# 147.92 564.9 m (2) jc s
# 175.14 47.242 m glx sp
# ...
# 727.19 564.9 m (X) jc s
# 744.57 47.242 m glx sp
# 751.22 564.9 m (Y) jc s
# 757.86 47.242 m glx sp
_num = pp.Combine(pp.Word(pp.nums) + "." + pp.Word(pp.nums)) | pp.Word(pp.nums)
_p_chr = (
    _num("START")
    + _num
    + pp.Literal("m")
    + pp.Literal("glx")
    + pp.Literal("sp")
    + _num
    + _num
    + pp.Literal("m")
    + pp.Literal("(")
    + (pp.Word(pp.nums) | pp.Literal("X") | pp.Literal("Y"))("CHR")
    + pp.Literal(")")
    + pp.Literal("jc")
    + pp.Literal("s")
)


def main():
    parser = argparse.ArgumentParser(
        description="Change the visual style of SMASH PostScript plots."
    )
    parser.add_argument(
        "--ps",
        "--input",
        type=argparse.FileType("r"),
        required=True,
        help="The PostScript file to read",
    )
    parser.add_argument(
        "--output",
        type=argparse.FileType("w"),
        required=True,
        help="The PostScript file to write",
    )
    parser.add_argument(
        "--highlight",
        type=str,
        metavar="CHR",
        choices=["".join(b) for b in product(_arms.keys(), ["p", "q"])],
        help="Highlight a given chromosome/arm. Examples: 1q, 7p, 8q",
    )
    parser.add_argument(
        "--color",
        type=str,
        metavar="COLOR",
        choices=_rgb.keys(),
        help="Color of highlight",
    )
    parser.add_argument("--title", type=str, help="New title for plots")

    args = parser.parse_args()

    # Get the chromosome start positions
    lines = args.ps.read()
    chr_start: Mapping[str, float] = {}
    starts = []
    for result, _, _ in _p_chr.scan_string(lines):
        start = float(result["START"])
        chr_start[result["CHR"]] = start
        starts.append(start)

    # Add chromosome mid and end positions
    # WARNING: This will always miscalculate the final chromosome, Y
    # This is because the true end is not parsed so we end up with a circle
    # chr_pos["Y"]["END"] == chr_pos["1"]["START"]
    chr_pos: Mapping[str, Mapping[str, float]] = {}
    for chrm, start in chr_start.items():
        try:
            end = starts[starts.index(start) + 1]
            mid = ((end - start) * _arms[chrm]["p"]) + start
            chr_pos[chrm] = {"START": start, "MID": mid, "END": end}
        except IndexError:
            pass

    start, end = 0, 0
    if args.highlight:
        arm = "p" if ("p" in args.highlight) else "q"
        chrm = args.highlight.replace(arm, "")

        if arm:
            if arm == "p":
                start = chr_pos[chrm]["START"]
                end = chr_pos[chrm]["MID"]
            elif arm == "q":
                start = chr_pos[chrm]["MID"]
                end = chr_pos[chrm]["END"]
        else:
            start = chr_pos[chrm]["START"]
            end = chr_pos[chrm]["END"]

    for line in lines.split("\n"):
        skip_next = False  # Default to not skipping

        if args.highlight:
            # Add a macro definition line with match
            if match := re.search("^\/p (.*) fp (.*)$", line):
                r, g, b = _rgb[args.color]
                ph = f"/ph {match.group(1)} {r / 255} {g / 255} {b / 255} setrgbcolor fp {match.group(2)}"
                args.output.write(f"{ph}\n")
                skip_next = False

            # Replace match with use of new macro
            if match := re.search("^(\d+\.?\d+) (\d+\.?\d+) p$", line):
                px = float(match.group(1))
                py = float(match.group(2))
                if (start <= px) and (px <= end):
                    args.output.write(f"{px} {py} ph\n")
                    skip_next = True

        if args.title:
            # Change the annotated title
            if re.search("^%%Title: .*$", line):
                args.output.write(f"%%Title: {args.title}\n")
                skip_next = True

            # Change the visual title
            if match := re.search("^(.*sf bk c \().*bins(\) jc s)\W*$", line):
                before = match.group(1)
                after = match.group(2)
                args.output.write(f"{before}{args.title}{after}\n")
                skip_next = True

        if not skip_next:
            args.output.write(f"{line}\n")


if __name__ == "__main__":
    main()
