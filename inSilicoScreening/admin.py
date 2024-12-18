from django.contrib import admin
from app.models import Target, Compound, KnownDrug

admin.site.register(Target)
admin.site.register(Compound)
admin.site.register(KnownDrug)