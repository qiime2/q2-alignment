import os.path
import shutil
import urllib.parse
from importlib.resources import files

import numpy as np
import q2templates
from pymsaviz import MsaViz, __version__
from q2_types.feature_data import AlignedDNAFASTAFormat
from skbio import DNA, TabularMSA

TEMPLATES = files("q2_alignment") / "_msa_visualizer" / "assets"


def _msa_stats(
    msa: TabularMSA,
) -> str:
    """Generates a table of statistics for the provided multiple sequence
    alignment in HTML.

    Args:
      msa: A TabularMSA object.

    Returns:
      A string containing the HTML tables with the statistics.
    """
    # Get the number of sequences and the total number of positions
    n_sequences, alignment_length = msa.shape

    # Calculate overall number of gaps
    ungapped_total = np.sum([len(seq.degap()) for seq in msa])
    accumulated_len = alignment_length * n_sequences
    gaps_total = _gaps(ungapped_total, accumulated_len)

    # Calculate overall GC content
    gc_per_sequence = [seq.gc_content() for seq in msa]
    overall_gc = np.mean(gc_per_sequence)

    html = f"""
<div class="row-mb-4">
    <table id="msa_stats", class="display">
        <thead>
            <tr>
                <th>Statistic</th>
                <th>Value</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td>Number of sequences</td>
                <td>{n_sequences}</td>
            </tr>
            <tr>
                <td>Total alignment length</td>
                <td>{alignment_length}</td>
            </tr>
            <tr>
                <td>Overall Gaps</td>
                <td>{gaps_total}</td>
            </tr>
            <tr>
                <td>Mean GC-content</td>
                <td>{overall_gc * 100.:.2f}%</td>
            </tr>
        </tbody>
    </table>
</div>
"""

    return html


def _gaps(ungapped_len: int, alignment_len: int):
    """Display the number of gaps per sequence in a formatted way.

    Args:
        ungapped_len: Ungapped-length of a sequence.
        alignment_len: The total length of a multiple sequence alignment or
          alignment.

    Returns:
        A string that shows the ungapped-length, the total alignment length,
        and the percentage of gaps in the format
        <gaps>/<total_len>(<percentage>%). For example:

        3/48(6%)
    """
    gaps = alignment_len - ungapped_len
    percentage = (gaps / alignment_len) * 100.0
    return f"{gaps}/{alignment_len}({percentage:.2f}%)"


def _per_sequence_stats(output_dir: str, msa: TabularMSA) -> str:
    """Generate a table of basic sequence statistics, BLAST URLs, and
    downloadable FASTA sequences for each sequences in the provided multiple
    sequence alignment in HTML.

    The BLAST URL that is constructed searches NCBI's nt database using a
    nucleotide query. One side effect of this function is that a FASTA file
    is written for each sequence contained within the provided MSA.

    Args:
      output_dir: Path to the directory to which the FASTA files will be
        written.
      msa: A TabularMSA object.

    Returns:
      A string containing the HTML table."""
    html_lines = [
        "<table id='per-sequence_stats', class='display'>",
        "<thead>",
        "<tr>",
        "<th>Sequence ID</th>",
        "<th>Ungapped-length</th>",
        "<th>Gaps</th>",
        "<th>GC-content</th>",
        "<th>Run</th>",
        "<th>Download as</th>",
        "</tr>",
        "</thead>",
        "<tbody>",
    ]
    # Get the total alignment length
    alignment_len = msa.shape[1]

    for seq in msa:
        seq_id = seq.metadata["id"]
        ungapped_seq = seq.degap()
        ungapped_len = len(ungapped_seq)
        gc_content = f"{seq.gc_content() * 100.:.2f}%"
        encoded_seq = urllib.parse.quote(str(ungapped_seq))

        # Construct BLAST URL
        base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        params_url = "?CMD=Put&PROGRAM=blastn&DATABASE=nt&QUERY="
        url = f"{base_url}{params_url}{encoded_seq}"

        # Write sequence to the FASTA file format
        fasta_basename = f"{seq_id}.fa"
        seq_fp = os.path.join(output_dir, fasta_basename)
        ungapped_seq.write(file=seq_fp, format="fasta")

        html_lines.append(
            f"<td>{seq_id}</td>"
            f"<td>{ungapped_len}</td>"
            f"<td>{_gaps(ungapped_len, alignment_len)}</td>"
            f"<td>{gc_content}</td>"
            f"<td><a href='{url}' target='_blank' class='btn btn-primary'>"
            f"BLAST</a></td>"
            f"<td><a href='{fasta_basename}' target='_blank' "
            f"class='btn btn-primary' "
            f"download='{fasta_basename}'>FASTA</a></td></tr>"
        )

    html_lines.append("</tbody>")
    html_lines.append("</table>")
    return "\n".join(html_lines)


def msa_visualizer(
    output_dir: str,
    alignment: AlignedDNAFASTAFormat,
    wrap_length: int = 80,
    show_count: bool = True,
    show_consensus: bool = False,
    dpi: int = 300,
) -> None:
    if alignment is None:
        raise ValueError(
            "Cannot visualize an empty multiple sequence alignment."
        )

    alignment_fp = str(alignment)

    msa = TabularMSA.read(alignment_fp, constructor=DNA)

    # Generate overall statistics for the multiple sequence alignment and
    # format into an HTML string
    msa_stats = _msa_stats(msa)

    # Generate a set of BLAST links for each sequence and format into an HTML
    # string
    sequence_stats = _per_sequence_stats(output_dir=output_dir, msa=msa)

    # Visualize multiple sequence alignment
    msa_visualization = MsaViz(
        msa=alignment_fp,
        format="fasta",
        wrap_length=wrap_length,
        show_count=show_count,
        show_consensus=show_consensus,
        color_scheme="Nucleotide",
    )

    # Save visualization to different file formats
    vis_filename = "msa_visualization"
    vis_basenames = []
    for extension in ["png", "jpg", "svg", "pdf"]:
        vis_basename = "{}.{}".format(vis_filename, extension)
        visualization_fp = os.path.join(output_dir, vis_basename)
        msa_visualization.savefig(os.path.join(visualization_fp), dpi=dpi)
        vis_basenames.append(vis_basename)

    png_file, jpg_file, svg_file, pdf_file = vis_basenames

    # Render the template
    index = os.path.join(TEMPLATES, "index.html")
    q2templates.render(
        index,
        output_dir,
        context={
            "msa_stats": msa_stats,
            "sequence_stats": sequence_stats,
            "png_file": png_file,
            "jpg_file": jpg_file,
            "svg_file": svg_file,
            "pdf_file": pdf_file,
            "wrap_length": wrap_length,
            "dpi": dpi,
            "pymsaviz_version": __version__,
        },
    )

    # Copy JavaScript and CSS files, required by DataTables, to the output
    # directory
    js = os.path.join(TEMPLATES, "js/dataTables.min.js")
    os.mkdir(os.path.join(output_dir, "js"))
    shutil.copy(js, os.path.join(output_dir, "dataTables.min.js"))

    css = os.path.join(TEMPLATES, "css/dataTables.min.css")
    os.mkdir(os.path.join(output_dir, "css"))
    shutil.copy(css, os.path.join(output_dir, "dataTables.min.css"))
