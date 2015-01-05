import sys

__all__ = ["estimate_snr", "mypolycos", "skytemp", "telescopes", "astro"]


def show_progress(iterator, width=0, tot=None, fmt="%d", show_number=False):
    """Wrap an iterator so that a progress counter is printed
        as we iterate.

        Inputs:
            iterator: The object to iterate over.
            width: The width of the progress bar.
                (Default: Don't show a progress bar, only the percentage)
            tot: The total number of iterations.
                (Default: Use len(iterator) to determine 
                    total number of iterations)

        Outputs:
            None
    """
    if tot is None:
        tot = len(iterator)
    old = -1
    curr = 1
    for toreturn in iterator:
        progfrac = curr/float(tot)
        progpcnt = int(100*progfrac)
        if progpcnt > old:
            neq = int(width*progfrac+0.5)
            nsp = width-neq
            bar = "["*bool(width) + \
                    "="*neq + \
                    " "*nsp + \
                    "]"*bool(width)
            old = progpcnt
            progpcntstr = fmt % progpcnt
            sys.stdout.write("     " + bar + " %s %% " % progpcntstr)
            if show_number:
                sys.stdout.write("(%d of %d)" % (curr, tot))
            sys.stdout.write("\r")
            sys.stdout.flush()
        curr += 1
        yield toreturn
    print "Done"


