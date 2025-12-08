from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),                 # main page
    path('run/', views.run_sampling, name='run_sampling'),  # ajax endpoint later
]
