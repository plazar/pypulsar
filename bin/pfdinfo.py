#!/usr/bin/env python
import argparse
import os.path

import prepfold


def main():
    for pfdfn in args.pfdfns:
        pfd = prepfold.pfd(pfdfn)
        lines = []
        if args.headers is not None:
            for header in args.headers:
                lines.append('# ' + args.sep.join(header.split(',')).decode('string-escape'))
        for attrs in args.attrs:
            vals = []
            for attr in attrs.split(','):
                if (attr[0] == '[') and (attr[-1] == ']'):
                    vals.append(attr[1:-1])
                else:
                    vals.append("%s" % getattr(pfd, attr))
            lines.append(args.sep.join(vals).decode('string-escape'))
        print "\n".join(lines)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get and format information " 
                                                 "from prepfold binary files.")
    parser.add_argument("pfdfns", nargs="+",
                        help="Prepfold binary files to grab information from.")
    parser.add_argument("-a", "--attr", dest='attrs', 
                        default=[], action='append',
                        help="Names of attributes to fetch. Multiple attributes "
                             "can be provided. Each should be separated by commas. "
                             "Multiple -a/--attrib flags can be specified. The "
                             "output from different -a/--attrib flags will be " 
                             "separated by newlines.")
    parser.add_argument("--sep", dest='sep', default=r'\t', 
                        help="Output separator for attributes on same line.")
    parser.add_argument("--header", dest='headers',
                        default=None, action='append',
                        help="Header text to show. Multiple header strings "
                             "can be provided. Each should be separated by commas. "
                             "Multiple --header flags can be specified, The "
                             "output from different --header flags will be separated "
                             "by newlines.")
    args = parser.parse_args()
    main()
