#!/usr/bin/env python
# -*- coding: utf8 -*-

__author__ = "Laurent Faucheux <faucheux@centre-cired.fr>"
__all__    = ['gwp_based_co2eq']


import numpy as np
import copy
import tools as ts

class gwp_based_co2eq(ts.cache):

	REFERENCES = {'[1]':
		      {'author(s)':'Solomon S., D. Qin, M. Manning, Z. Chen, M. Marquis, K.B. Averyt, M. Tignor, H.L. Miller',
		       'year':2007,
                       'publisher':'Cambridge University Press',
		       'review':None,
		       'volume':None,
		       'numero':None,
		       'issue':None,
		       'pages':'211-213',
		       'url':'http://www.ipcc.ch/pdf/assessment-report/ar4/wg1/ar4_wg1_full_report.pdf',
		       'title':'Climate Change 2007: The Physical Science Basis. Contribution of Working Group I to the Fourth Assessment Report of the Intergovernmental Panel on Climate Change',
		       },
                      '[2]':
		      {'author(s)':u'Annie Levasseur, Pascal Lesage, Manuele Margni, Louise Deschênes, Réjean Samson',
		       'year':2010,
                       'publisher':'ACS Publications',
		       'review':'Environmental Science & Technology',
		       'volume':'44',
		       'numero':None,
		       'issue':'8',
		       'pages':'3169–3174',
		       'url':'http://pubs.acs.org/doi/abs/10.1021/es9030003',
		       'title':'Considering Time in LCA: Dynamic LCA and Its Application to Global Warming Impact Assessments',
		       },
		      }

	GHGS_BASE_DATA = {
				'CO2'    : {'a'           :{0:0.217, 1:0.259, 2:0.338, 3:0.186},
					    'tau'         :{1:172.9, 2:18.51, 3:1.186},
					    'GWP100'      :1.,
					    'AGWP100'     :8.69e-14,
					    'trajectories':{'AGWP':{},'GWP':{}},
					    },
				'N2O'    : {'tau'         :{'i':114.},
					    'GWP100'      :298.,
					    'trajectories':{'AGWP':{},'GWP':{}},
					    },
				'CH4'    : {'tau'         :{'i':12.},
					    'GWP100'      :25.,
					    'trajectories':{'AGWP':{},'GWP':{}},
					    },
				}
	OTHER_GHGS = ['N2O','CH4']


	def __init__(self,first_year,project_horizon,GWP_horizon=100,static=False):
		u"""
		Computation method which can consider time in life cycle analysis,
		based on the IPCC WG1's AR4 (2007) [1]_ and Levasseur et al (2010) [2]_.

		Parameters
		----------
		first_year      : int
		project_horizon : int
		GWP_horizon     : int, optional (100 by default)
		static          : bol, optional (False by default)
				  if True, conversions coefficients from weight of GHG to co2eq weight are
				  time-constant.

		Public attributes
		-----------------
		REFERENCES     : dict
				 applied background theory
		GHGS_BASE_DATA : dict
				 background data used for calculations
		OTHER_GHGS     : list
				 Other than CO2 greenhouse gases
		ghgs_data      : dict
				 GWP_horizon specific data used for calculations

		References
		----------

		.. [1] Solomon S., D. Qin, M. Manning, Z. Chen, M. Marquis, K.B. Averyt, M. Tignor, H.L. Miller (2007)
		    "Climate Change 2007: The Physical Science Basis. Contribution of Working Group I to the Fourth Assessment
		    Report of the Intergovernmental Panel on Climate Change".
		    Cambridge University Press 2007, Chapter 2, Section 2.10.2, Table 2.14, pp 211-213.

		.. [2] Annie Levasseur, Pascal Lesage, Manuele Margni, Louise Deschênes, Réjean Samson (2010)
		    "Considering Time in LCA: Dynamic LCA and Its Application to Global Warming Impact Assessments".
		    Environmental Science & Technology, 2010, 44 (8), pp 3169–3174.

		Example
		-------
		>>> dyn_gwp20  = gwp_based_co2eq(first_year=2020,
						 project_horizon=5,
						 GWP_horizon=20,
						 static=False)
		>>> print(dyn_gwp20.OTHER_GHGS)
		['N2O','CH4']
		"""

		super(gwp_based_co2eq,self).__init__()
		self.first_year       = first_year
		self.project_horizon  = project_horizon
		self.last_year        = first_year + self.project_horizon
		self.GWP_horizon      = float(GWP_horizon) if GWP_horizon else float(self.project_horizon)
		self.static           = static

	@ts.cache._property
	def ghgs_data(self):

		_bd = copy.deepcopy(self.GHGS_BASE_DATA)

		# < AGWPs >
		for ghg in self.OTHER_GHGS:
			_bd[ghg]['AGWP100'] = _bd['CO2']['AGWP100']\
					      * _bd[ghg]['GWP100']

		_bd['CO2']['integral100'] = _bd['CO2']['a'][0]*(1e2-1)\
					    + sum(
						    map(lambda i:_bd['CO2']['a'][i]*_bd['CO2']['tau'][i]
							*(1-np.exp((-1e2)/_bd['CO2']['tau'][i])),range(1,3+1)
							)
						    )
		# < INTEGRALS >
		for ghg in self.OTHER_GHGS:
			_bd[ghg]['integral100'] = _bd[ghg]['tau']['i']\
						  *(1-np.exp((-1e2)/_bd[ghg]['tau']['i']))
		# < ais >
		for ghg in _bd.keys():
			_bd[ghg]['ai100'] = _bd[ghg]['AGWP100']/_bd[ghg]['integral100']

		# < AGWPi(t) (W/m2/kg) >
		for _t_ in range(0,self.project_horizon)[:(1 if self.static else None)]:
			year = self.first_year + _t_
			_bd['CO2']['trajectories']['AGWP'][year] = _bd['CO2']['ai100']\
								   * (_bd['CO2']['a'][0]*(self.GWP_horizon-_t_-1)\
								      + sum(
									      map(lambda i:_bd['CO2']['a'][i]\
										  *_bd['CO2']['tau'][i]\
										  *(1-np.exp((_t_-self.GWP_horizon)/_bd['CO2']['tau'][i])),\
										  range(1,3+1)
										  )
									      )
								      )
			for ghg in self.OTHER_GHGS:
				_bd[ghg]['trajectories']['AGWP'][year] = _bd[ghg]['ai100']\
									 * _bd[ghg]['tau']['i']\
									 *(1-np.exp(((_t_-self.GWP_horizon)/_bd[ghg]['tau']['i'])))
			for ghg in _bd.keys():
				_bd[ghg]['trajectories']['GWP'][year] = _bd[ghg]['trajectories']['AGWP'][year]\
									/_bd['CO2']['trajectories']['AGWP'][self.first_year]
		return _bd

	def co2eq_yields_trajectory_computer(self,ghgs_yield,as_row_array=True):
		"""
		Core method which computes the total quantity of CO2eq per weight unit of emitting source.

		Parameters
		----------
		ghgs_yield : dict
		  Yields related to weight of each emitted ghg per weight of emitting source

		as_row_array : bol, optional (True by default)

		Example 1
		---------
		>>> dyn_gwp20  = ...
		>>>
		>>> ghgs_weight_per_weight_of_output_inventory_flow = {'CO2':1.,'N2O':.0,'CH4':.0}
		>>>
		>>> co2eq_traj = dyn_gwp20.co2eq_yields_trajectory_computer(ghgs_weight_per_weight_of_output_inventory_flow,as_row_array=False)
		>>> print(co2eq_traj['as_array'])
		[[ 1.        ]
		 [ 0.95764081]
		 [ 0.91469171]
		 [ 0.87112496]
		 [ 0.82691128]]
		>>> print(co2eq_traj['as_dict'])
		{2024: 0.82691127746144444, 2020: 1.0, 2021: 0.95764080833063492, 2022: 0.91469171438570718, 2023: 0.87112496115582216}


		Example 2
		---------
		>>> dyn_gwp20  = ...
		>>>
		>>> co2eq_traj = dyn_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':1.,'CH4':.0})
		>>> print(co2eq_traj['as_array'])                                           ^^
		[[ 292.33637282  278.90543843  265.35617058  251.68752668  237.89845498]]


		Example 3
		---------
		>>> dyn_gwp20  = ...
		>>>
		>>> co2eq_traj = dyn_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':.0,'CH4':1.})
		>>> print(co2eq_traj['as_array'])                                                    ^^
		[[ 72.2209832   70.75950679  69.17102216  67.44449179  65.56791893]]


		Example 4
		---------
		>>> sta_gwp20  = gwp_based_co2eq(first_year=2020,
						 project_horizon=5,
						 GWP_horizon=20,
						 static=True)
		>>>                                     ^^^^
		>>> co2eq_traj = sta_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':.0,'CH4':1.})
		>>> print(co2eq_traj['as_array'])
		[[ 72.2209832  72.2209832  72.2209832  72.2209832  72.2209832]]

		"""
		trajectory = {}
		for _t_ in range(self.project_horizon):
			year = self.first_year + _t_
			# < GWPi(t) (kgCO2eq/kgi) & CO2eq(t) >
			trajectory[year] = 0.
			for ghg in self.ghgs_data.keys():
				trajectory[year] += self.ghgs_data[ghg]['trajectories']['GWP'][self.first_year if self.static else year]\
						    * ghgs_yield[ghg]
		return {'as_dict':trajectory,
			'as_array':np.array(
					 map( lambda item:item[1], sorted(trajectory.items(),key=lambda item:item[0])
					     )
					    ).reshape((1 if as_row_array else self.project_horizon,
						       self.project_horizon if as_row_array else 1))}


if __name__ == '__main__':

	cls = lambda H,S : gwp_based_co2eq(first_year=2020,project_horizon=5,GWP_horizon=H,static=S)

	print("\nExample 1, 2 & 3\n"+16*"-")

	dyn_gwp20  = cls(20,False)
	print(dyn_gwp20.co2eq_yields_trajectory_computer({'CO2':1.,'N2O':.0,'CH4':.0}).items())
	print(dyn_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':1.,'CH4':.0})['as_array'])
	print(dyn_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':.0,'CH4':1.})['as_array'])


	print("\nExample 4\n"+9*"-")

	sta_gwp20  = cls(20,True)
	print(sta_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':.0,'CH4':1.})['as_array'])
