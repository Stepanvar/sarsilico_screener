from django.urls import path
from app import views

urlpatterns = [
    path('', views.home_view, name='home'),
    path('about/', views.about_view, name='about'),
    path('contact/', views.contact_view, name='contact'),
    path('login/', views.login_view, name='login'),
    path('logout/', views.logout_view, name='logout'),
    path('smiles-input/', views.smiles_input_view, name='smiles_input'),
    path('pdb-upload/', views.pdb_upload_view, name='pdb_upload'),
    path('docking/', views.docking_view, name='docking'),
    path('docking/results/<uuid:job_id>/', views.docking_results_view, name='docking_results'),
    path('similarity/', views.similarity_analysis_view, name='similarity_analysis'),
    path('similarity/results/<uuid:job_id>/', views.similarity_results_view, name='similarity_results'),
    path('sequence-conversion/', views.sequence_conversion_view, name='sequence_conversion'),
    path('select/target/', views.target_selection_view, name='target_selection'),
    path('register/', views.register_view, name='register'),
]
