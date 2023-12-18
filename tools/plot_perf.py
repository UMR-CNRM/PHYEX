#!/usr/bin/env python3

"""
This script plots the data contained in the performance files obtained with the --perf option
"""

import matplotlib.pyplot as plt
import numpy
import pandas

class Perf():
    def __init__(self, perffile):
        """
        :param perffile: text file with each line having the form "commit model case time"
        """
        self._df = df = pandas.read_csv(perffile, sep=' ',
                                        names=['commit', 'model', 'case', 'time'])
    def plotPerf(self, outfile, model=None, title=None, num=None):
        """
        :param outfile: output file
        :param model: None to plot each model on a subplot or the model to plot
        :param title: custom title to use (%M will be replaced by the model name)
        :param num: plot only last num values (None to plot all values)
        """
        models = [model] if model is not None else sorted(set(self._df['model']))
        fig, ax = plt.subplots(nrows=len(models), sharex=True, sharey=True, figsize=(8, 8 * len(models)))
        if len(models) == 1:
            ax = [ax]
    
        #Ordered commit list common to all models
        commits = []
        for commit in self._df['commit']:
            if commit not in commits:
                commits.append(commit)
        if num is not None:
            commits = commits[-num:]
    
        df = self._df.groupby('model')
        for igrpM, grpM in enumerate(models):
            if title is None:
                if len(models) == 1:
                    ax[igrpM].set_title('Mean elapsed computational time')
                else:
                    ax[igrpM].set_title('Mean elapsed computational time for ' + grpM)
            else:
                ax[igrpM].set_title(title.replace('%M', grpM))
            ax[igrpM].set_ylabel('time (ms/gp)')
            ax[igrpM].set_yscale('log')
    
            dfp = df.get_group(grpM).groupby('case')
            for grp in dfp.groups:
                #Build time serie with possible missing value and mean aggregation if needed
                time = []
                for commit in commits:
                    f = dfp.get_group(grp)['commit'] == commit
                    #discard negative values
                    l = [numpy.nan if t < 0. else t for t in dfp.get_group(grp)[f]['time']]
                    time.append(numpy.nan if len(l) == 0 else numpy.ma.array(l).mean())
                ax[igrpM].plot(range(len(commits)), numpy.ma.array(time), 'o-', label=grp)
            if igrpM == len(df.groups) - 1:
                ax[igrpM].set_xlabel('PHYEX version')
                ax[igrpM].set_xticks(range(len(commits)))
                ax[igrpM].set_xticklabels(commits, rotation=45, ha='right')
            ax[igrpM].legend()
        fig.tight_layout()
        fig.savefig(outfile)

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
    args = parser.parse_args()
    perf = Perf(args.PERF_FILE)
    if args.plot is not None:
        perf.plotPerf(args.plot, args.model, args.title, args.num)
    if args.listModels:
        print(' '.join(perf.listModels()))
        
                        
