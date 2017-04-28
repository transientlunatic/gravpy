"""
Provides data from the ATNF catalogue.
"""
import urllib
import numpy as np
import os
import pandas as pd
import datetime
import astropy.table
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

def get_atnf():
        """
        Contact the ATNF pulsar catalogue, returning an array of data.
        Each row corresponding to one pulsar, with columns in the format:
        """

        tnow = datetime.datetime.today()
        try:
            tdata = datetime.datetime.fromtimestamp(
                    os.path.getmtime('pulsar_data.txt')
                    )
        except:
            tdata = datetime.datetime(1977,1,1)
        dt = tnow - tdata
        if dt.days > 20:

            try:
                #URL to get |NAME|PSRJ|RAJ|DECJ|POSEPOCH|F0|F1|F2|PEPOCH|DM|W50|W10|S400|S1400|SPINDX|PSRTYPE|NGLT|
                url2 = 'http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.54&Name=Name&JName=JName&RaJ=RaJ&DecJ=DecJ&PosEpoch=PosEpoch&F0=F0&F1=F1&F2=F2&PEpoch=PEpoch&DM=DM&DIST=DIST&W50=W50&W10=W10&S400=S400&S1400=S1400&SPINDX=SPINDX&Type=Type&NGlt=NGlt&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj/decj&radius=&coords_1=&coords_2=&style=Short+without+errors&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=40&table_bottom.y=24'
                Hurl2='#NAME|PSRJ|RAJ|DECJ|POSEPOCH|F0|F1|F2|PEPOCH|DM|W50|W10|S400|S1400|SPINDX|PSRTYPE|NGLT\n'
                sock = urllib.urlopen(url2)
                data = sock.read()
                sock.close()


                data = data.split('<pre>')[1]
                data = data.split('</pre>')[0]
                data = data.splitlines()[5:-1]

                header = "Name|PSRJ|RAJ|DECJ|POSEPOCH|F0|F1|F2|PEPOCH|DM|DIST|W50|W10|S400|S1400|SPINDX|PSRTYPE|NGLT".split('|')

                df = pd.DataFrame(columns=header)

                for b in data:
                    b = b.split()
                    try:
                        df = df.append(pd.Series(b[1:], index=header, name=b[0]))
                    except:
                        pass
                del(df[''])
            except:
                pass
        
            df.to_csv('pulsar_data.txt', header=header)
        else:
            df = pd.read_csv('pulsar_data.txt', na_values='*')
        
        a = astropy.table.Table.from_pandas(df)
        
        a['RAJ'].unit = u.hourangle
        a['DECJ'].unit = u.deg
        a['F0'].unit = u.hertz
        a['F1'].unit = u.hertz / u.second
        a['F2'].unit = u.hertz / u.second**2
        a['DIST'].unit = u.kiloparsec
        a['W50'].unit = u.millisecond
        a['W10'].unit = u.millisecond
        a['DM'].unit = u.centimeter**(-3) * u.parsec
        #a['EDOT'].unit = u.erg / u.second
        a.add_index('Unnamed: 0')
        a.rename_column('Unnamed: 0', '#')
        a.add_index('PSRJ')
        a['POS'] = SkyCoord(a['RAJ'], a['DECJ'], unit=(a['RAJ'].unit, a['DECJ'].unit))
        return a