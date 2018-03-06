#!/usr/bin/env python
# -*- coding: utf8 -*-

__author__ = "Laurent Faucheux <faucheux@centre-cired.fr>"
__all__    = ['cache']


class cache(object):

        def __init__(self,*args,**kwargs):
                """ Homemade cache class which aims at being inherited """
                self._cache = {}

	@classmethod
	def _property(cls,meth):
		@property
		def __property(cls,*args,**kwargs):
			meth_name = meth.__name__
			if meth_name not in cls._cache:
				cls._cache[meth_name] = meth(cls,*args,**kwargs)
			return cls._cache[meth_name]
		return __property
