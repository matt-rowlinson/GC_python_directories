#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

def read_nitrate_dataset(convert=False):
    path='/mnt/lustre/users/mjr583/NCAS_CVAO/CVAO_datasets/'
    df = pd.read_csv(path+'fomba_CV_aerosol_2014.csv', skiprows=37, index_col=1)
    df.index = pd.to_datetime(df.index, format="%Y-%m-%dT%H:%M")
    df.index = df.index.tz_localize('UTC').tz_convert('Atlantic/Cape_Verde')

    df = df[ [df.columns[7]]  ]
    df.columns = ['NO3 / ug/m3']
    if convert:

        #df = df / ((28.9644 / 62.0049) )

        df = 24.45 * df / 62.0049
        df.columns = ['NO3 / ppb']

    return df

def plot_nitrate():
    plot = df.plot()
    fig = plot.get_figure()
    fig.savefig("output.png")
    return

def main():
    df = read_nitrate_dataset(convert=False)
    print( df )
    df = read_nitrate_dataset(convert=True)
    print( df )



if __name__ == "__main__":
    main()

