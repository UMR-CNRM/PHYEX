#!/usr/bin/env python3

"""
This script plots the data contained in the performance files obtained with the --perf option
"""

import re
import matplotlib.pyplot as plt
import numpy
import pandas

class Perf():
    """
    This class manipulates performance files
    """
    def __init__(self, perffile):
        """
        :param perffile: text file with each line having the form "commit model case time"
        """
        self._df = pandas.read_csv(perffile, sep=' ',
                                   names=['commit', 'model', 'case', 'time', 'time2'])

    def plotPerf(self, outfile, model=None, title=None, num=None, allTimes=False):
        """
        :param outfile: output file
        :param model: None to plot each model on a subplot or the model to plot
        :param title: custom title to use (%M will be replaced by the model name)
        :param num: plot only last num values (None to plot all values)
        :param allTimes: True to plot alternative metrics
        """
        models = [model] if model is not None else sorted(set(self._df['model']))
        fig, ax = plt.subplots(nrows=len(models), sharex=True, sharey=True,
                               figsize=(8, 8 * len(models)))
        if len(models) == 1:
            ax = [ax]
 
        #Ordered commit list common to all models
        commits = []
        for commit in self._df['commit']:
            if commit not in commits:
                commits.append(commit)
        if num is not None:
            commits = commits[-num:]
        shortCommits = [self.shortenCommit(c) for c in commits]

        df = self._df.groupby('model')
        for igrpM, grpM in enumerate(models):
            if title is None:
                if len(models) == 1:
                    ax[igrpM].set_title('Mean elapsed computational time')
                else:
                    ax[igrpM].set_title('Mean elapsed computational time for ' + grpM)
            else:
                ax[igrpM].set_title(title.replace('%M', grpM))
            ax[igrpM].set_ylabel('time')
            ax[igrpM].set_yscale('log')

            dfp = df.get_group(grpM).groupby('case')
            for grp in dfp.groups:
                #Build time serie with possible missing value and mean aggregation if needed
                time = []
                time2 = []
                for commit in commits:
                    f = dfp.get_group(grp)['commit'] == commit
                    #discard negative values
                    l = [numpy.nan if t < 0. else t for t in dfp.get_group(grp)[f]['time']]
                    time.append(numpy.nan if len(l) == 0 else numpy.ma.array(l).mean())
                    l = [numpy.nan if t < 0. else t for t in dfp.get_group(grp)[f]['time2']]
                    time2.append(numpy.nan if len(l) == 0 else numpy.ma.array(l).mean())
                p = ax[igrpM].plot(range(len(commits)), numpy.ma.array(time), 'o-', label=grp)
                if allTimes:
                    ax[igrpM].plot(range(len(commits)), numpy.ma.array(time2), 'o:', color=p[0].get_color())
            if igrpM == len(models) - 1:
                ax[igrpM].set_xlabel('PHYEX version')
                ax[igrpM].set_xticks(range(len(commits)))
                ax[igrpM].set_xticklabels(shortCommits, rotation=45, ha='right')
            ax[igrpM].legend()
        fig.tight_layout()
        fig.savefig(outfile)

    @staticmethod
    def shortenCommit(commit):
        """
        :param commit: full commit SHA
        :return: shorten version of commit SHA
        """
        return commit[:7] if re.match(r'^[0-9a-f]{40}$', commit) else commit

    def listModels(self):
        """
        :result: list of models present in the file
        """
        return sorted(set(self._df['model']))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Plot performance file')
    parser.add_argument('PERF_FILE', type=str,
                        help='file containing the performance statistics')
    parser.add_argument('--plot', type=str, default=None,
                        help='output plot file')
    parser.add_argument('--model', type=str, default=None,
                        help='plot only model MODEL')
    parser.add_argument('--title', type=str, default=None,
                        help="Plot title, %%M will be replaced by model name")
    parser.add_argument('--num', metavar='N', type=int, default=None,
                        help="Plot only last N values")
    parser.add_argument('--listModels', default=False, action='store_true',
                        help="returns the list of models present in the performance file")
    parser.add_argument('--allTimes', default=False, action='store_true',
                          help="to also plot the alternative metrics")
    args = parser.parse_args()
    perf = Perf(args.PERF_FILE)
    if args.plot is not None:
        perf.plotPerf(args.plot, args.model, args.title, args.num, args.allTimes)
    if args.listModels:
        print(' '.join(perf.listModels()))
