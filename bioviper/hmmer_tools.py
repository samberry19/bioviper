from __init__ import *
import subprocess
import Bio

def hmmalign(sequences, hmm, outfilename='tmp.sto', reload=True):

    if type(sequences)==str:
        seqfilename = sequences
    else:
        Bio.SeqIO.write(sequences, "tmpseqs.fasta", "fasta")
        seqfilename = "tmpseqs.fasta"

    if type(hmm)==str:
        hmmfilename = hmm
    elif isinstance(hmm, MultipleSequenceAlignment):
        msa = hmm
        hmmbuild(msa, "tmp.hmm")
        hmmfilename = "tmp.hmm"

    cmd = "hmmalign -o "+outfilename+" " + hmmfilename + " " + seqfilename
    subprocess.call("module load hmmer;"+cmd, shell=True)

    if reload==True:
        ali = readAlignment(outfilename, "stockholm")
        return ali

def hmmbuild(msa, hmm_filename, aliname="tmp.fa"):

    msa.save(aliname)
    subprocess.call("module load hmmer; hmmbuild " + hmm_filename + " "  + aliname, shell=True)

#def reformat(input_filename,

    # Reformat an alignment using esl-reformat without explicitly loading it into memory -
        # super fast

default_hmmer_bin = "/n/gaudet_lab/programs/hmmer-bin/"

def esl_alimask(alignment, save_as='tmp.sto', output = "Default", hmmer_dir=default_hmmer_bin,
               p = True, ppcons = True, ppcons_thresh = 0.5, g = False, gapthresh = 0.5,
               final_format = "stockholm", load = True):

    """ This serves as a wrapper function to call the esl tool alimask, which is great at very quickly subsetting columns
        of a huge alignment based on gap frequencies or posterior probability values. While these are all possible with the
        MultipleSequenceAlignment object, doing so for alignments with many columns is slow and requires a huge amount of
        memory, while this can do so in seconds.

        By default, this function masks based on a 50% consensus probability threshold, but can support a variety of other
        functionalities:
            p = True/False (Default: TRUE) controls whether a probability threshold is used
            ppcons = True/Fase (Default: TRUE) controls whether the probability threshold is based on the consensus probability
                annotation for each column.
            ppcons_thresh = <float> (0-1) (Default: 0.5) sets the consensus probability threshold.

            g = True/False (Default: FALSE): masks by fraction of gaps, either in addition to or instead of the probability threshold.
            gapthresh = <float> (Default: 0.5) sets the gap fraction, between 0 and 1

        """

    if type(alignment)==str:
        alifilename = sequences

    else:
        alignment.save(save_as, "stockholm")
        alifilename = save_as

    if output == "Default":
        output = ''.join(alifilename.split('.')[:-1]) + '_masked.sto'

    cmd = hmmer_dir + os.sep + "esl-alimask -o " + output

    if p:
        cmd = cmd + ' -p'

        if ppcons:
            cmd = cmd + ' --ppcons ' + str(ppcons_thresh)

    if g:
        cmd = cmd + ' -g'

        if gap_thresh != 0.5:
            cmd = cmd + ' --gapthresh ' + str(gapthresh)

    subprocess.call(cmd, shell=True)

    if final_format != 'stockholm':

        ali = reformat(output, "afa")

        if load:
            return ali

    elif load:
        return readAlignment(output, "stockholm")
