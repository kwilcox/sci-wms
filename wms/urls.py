'''
COPYRIGHT 2010 RPS ASA

This file is part of SCI-WMS.

    SCI-WMS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SCI-WMS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SCI-WMS.  If not, see <http://www.gnu.org/licenses/>.
'''

from django.conf.urls import patterns, url

urlpatterns = patterns( '',

                        url(r'^index', 'wms.views.index'),
                        url(r'^$', 'wms.views.index'),

                        url(r'^documentation/', 'wms.views.documentation', name='documentation'),

                        # Datasets
                        url(r'^datasets$',  'wms.views.datasets'),
                        url(r'^datasets/$', 'wms.views.datasets'),
                        url(r'^datasets/(?P<dataset>.*)/update', 'wms.views.update_dataset', name="update_dataset"),
                        url(r'^datasets/(?P<dataset>.*)/', 'wms.views.wms', name="dataset"),

                        url(r'^standard_names', 'sciwms.apps.wms.views.standard_names'),

                        # Colormaps
                        url(r'^colormaps', 'sciwms.apps.wms.views.colormaps'),

                        # Clients
                        url(r'^openlayers/(?P<filepath>.*)', 'wms.views.openlayers'),
                        url(r'^simple', 'wms.views.simpleclient', name='simpleclient'),
                        url(r'^leaflet', 'wms.views.leafletclient', name='leafletclient'),

                        url(r'^add_dataset', 'wms.views.add'),  # This is a POST based view
                        url(r'^add_to_group', 'wms.views.add_to_group'),
                        url(r'^remove_dataset', 'wms.views.remove'),
                        url(r'^remove_from_group', 'wms.views.remove_from_group'),

                        url(r'^groups/(?P<group>.*)/wmstest/$', 'wms.views.grouptest'),
                        url(r'^groups/(?P<group>.*)/wmstest', 'wms.views.grouptest'),
                        url(r'^groups/(?P<group>.*)/', 'wms.views.groups'),
                        url(r'^groups/(?P<group>.*)', 'wms.views.groups')
                    )
