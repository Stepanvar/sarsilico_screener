from django.contrib import admin
from app.models import PredefinedProteinTarget, Compound, ScreeningJob

admin.site.register(PredefinedProteinTarget)
admin.site.register(Compound)
admin.site.register(ScreeningJob)