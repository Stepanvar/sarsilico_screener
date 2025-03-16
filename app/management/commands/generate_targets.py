# app/management/commands/generate_targets.py
import os
from pathlib import Path
from django.core.management.base import BaseCommand
from django.conf import settings
from app.models import PredefinedProteinTarget

class Command(BaseCommand):
    help = 'Delete existing targets and generate new ones from PDB files'

    def handle(self, *args, **options):
        # Delete existing targets
        PredefinedProteinTarget.objects.all().delete()
        self.stdout.write(self.style.SUCCESS('Deleted all existing targets'))

        # Path to PDB files
        pdb_dir = Path('/home/s-zuev/inSilicoScreening/media/pdb_files')
        
        # Counter for created targets
        created_count = 0

        try:
            # Iterate through all PDB files
            for file in pdb_dir.glob('*.pdb'):
                # Extract base name without extension
                name = file.stem
                
                if not name:
                    self.stdout.write(self.style.WARNING(
                        f'Skipping invalid filename: {file.name}'
                    ))
                    continue

                # Create new target
                PredefinedProteinTarget.objects.create(
                    name=name,
                    description=f"Auto-generated from PDB file {file.name}",
                    pdb_file=f'pdb_files/{file.name}'
                )
                created_count += 1

            self.stdout.write(self.style.SUCCESS(
                f'Successfully created {created_count} targets'
            ))

        except FileNotFoundError:
            self.stdout.write(self.style.ERROR(
                f'PDB directory not found: {pdb_dir}'
            ))
        except Exception as e:
            self.stdout.write(self.style.ERROR(
                f'Error processing files: {str(e)}'
            ))
