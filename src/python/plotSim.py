import argparse
import os
import afinesanalysis.afinesanalysis as aa
import afinesanalysis.plotting as pltaf

parser = argparse.ArgumentParser()

parser.add_argument('directory', help='directory of simulation to be plotted')
parser.add_argument('-viz', '--visualize', type=str, choices=['all', 'filaments'],
                    help='What to plot. Options are filaments or all, defaults to \
                    filaments', default='filaments')
parser.add_argument('-dt', help='Plot every dt-th frame. Defaults to 10',
                    type=int, default=10)
parser.add_argument('-df', '--dfilament', help='Plot every dfilament-th filament. \
                    Defaults to 1, i.e. plot every filament', type=int, default=1)
parser.add_argument('-dp', '--dpmotor', help='Plot every dpmotor-th pmotor. \
                    Defaults to 1', type=int, default=1)
parser.add_argument('-da', '--damotor', help='Plot every damotor-th amotor. \
                    Defaults to 1', type=int, default=1)

args = parser.parse_args()
# print(args)
os.chdir(args.directory)

configs = aa.readConfigs(os.path.join(args.directory, 'data', 'config_full.cfg'))
filamentData = aa.readData(os.path.join(args.directory, 'txt_stack', 'actins.txt'),
                           configs)

if args.visualize == 'filaments':
    pltaf.filaments(filamentData, configs, dt=args.dt, dfilament=args.dfilament,
                    savepath=os.path.join(args.directory, 'filaments'))
elif args.visualize == 'all':
    amotorData = aa.readData(os.path.join(args.directory, 'txt_stack', 'amotors.txt'),
                             configs)
    pmotorData = aa.readData(os.path.join(args.directory, 'txt_stack', 'pmotors.txt'),
                             configs)
    pltaf.all(filamentData, pmotorData, amotorData, configs, args.dt, args.dfilament,
              args.dpmotor, args.damotor,
              savepath=os.path.join(args.directory, 'all'))
