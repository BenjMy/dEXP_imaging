# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 12:25:29 2020

@author: Benjamin
"""

import pandas as pd

df = pd.DataFrame(np.random.rand(5, 6), columns=list('abcdef'))


df['a'].split([4],expand=True)

df.a.str.split(expand=True,)
