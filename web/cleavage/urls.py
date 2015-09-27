from django.conf.urls import patterns, url
from django.views.generic import TemplateView

import views

urlpatterns = patterns('',
    
    url(r'^$',
        TemplateView.as_view(template_name = 'home.html'),
    ),
    
    url(r'^cleavage-prediction/$',
        views.cleavage_prediction,
    ),
)
