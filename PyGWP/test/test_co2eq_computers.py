#!/usr/bin/env python
# -*- coding: utf8 -*-

__author__ = "Laurent Faucheux <faucheux@centre-cired.fr>"
__all__    = ['test_gwp_based_co2eq']

import unittest
import numpy as np
from PyGWP import gwp_based_co2eq


class test_gwp_based_co2eq(unittest.TestCase):
	_tested_inst = gwp_based_co2eq(first_year=0,
				       project_horizon=1,
				       GWP_horizon=100,
				       static=True)

	def test_mandatory_REFERENCES_keys(self):
		"""
		Method which tests that non-facultative references keys are not None
		"""
		Trues = []
		for mandatory_key in ['author(s)','year','title']:
			Trues += map(lambda d:d.get(mandatory_key),
				      self._tested_inst.REFERENCES.values()
				      )
		self.assertTrue(all(Trues))

	def test_GHGS_BASE_DATA_structure(self):
		"""
		Method which tests that each GHGS_BASE_DATA
		greenhouse gase dict contains a GWP100 key
		"""
		self.assertTrue(all(map(lambda d:d.has_key('GWP100'),
					self._tested_inst.GHGS_BASE_DATA.values()
					)
				    )
				)

	def test_N2O_RGWP(self):
		""" Method which tests the N2O recomputed relative GWP100 """
		array_res = self._tested_inst\
			    .co2eq_yields_trajectory_computer({'CO2':.0,
							       'N2O':1.,
							       'CH4':.0}\
							      )['as_array']
		array_ref = np.array(
			[[ self._tested_inst.GHGS_BASE_DATA['N2O']['GWP100'] ]]
			)
		np.testing.assert_equal(np.round(array_res,12),
					np.round(array_ref,12))

	def test_CH4_RGWP(self):
		""" Method which tests the CH4 recomputed relative GWP100 """
		array_res = self._tested_inst\
			    .co2eq_yields_trajectory_computer({'CO2':.0,
							       'N2O':.0,
							       'CH4':1.}\
							      )['as_array']
		array_ref = np.array(
                        [[ self._tested_inst.GHGS_BASE_DATA['CH4']['GWP100'] ]]
                        )
		np.testing.assert_equal(np.round(array_res,12),
					np.round(array_ref,12))

if __name__ == '__main__': unittest.main()
