# app/services/similarity_checker.py
"""
Drug repurposing checks implementing insilico.txt section 3 requirements
Handles Morgan fingerprints, TTD API queries, and CoviDrug scraping
"""

import logging
import requests
from rdkit import Chem
from rdkit.Chem import DataStructs, rdMolDescriptors
from django.conf import settings
from typing import Dict, Optional
from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)

class SimilarityChecker:
    """
    Drug repurposing validation system implementing:
    - Morgan fingerprint generation (2048 bits)
    - TTD API similarity search
    - CoviDrug web scraping
    """
    
    def __init__(self):
        self.ttd_endpoint = settings.TTD_API_ENDPOINT
        self.covidrug_base = settings.COVIDRUG_SCRAPE_URL
        self.similarity_threshold = settings.SIMILARITY_THRESHOLD

    def check_similarity(self, smiles: str) -> Dict:
        """
        Main workflow for drug repurposing check
        Returns: Dict with similarity data and warnings
        """
        try:
            # Generate validated fingerprint
            fp = self.generate_fingerprint(smiles)
            if not fp:
                return {'error': 'Invalid compound for fingerprint generation'}

            # Query Therapeutic Target Database
            ttd_results = self.query_ttd(fp)
            if 'error' in ttd_results:
                return ttd_results

            # Cross-reference with CoviDrug if threshold met
            if ttd_results['max_similarity'] >= self.similarity_threshold:
                covi_data = self.scrape_covdrug(ttd_results['best_match_id'])
                return {
                    **ttd_results,
                    'covi_data': covi_data,
                    'warning': True,
                    'message': 'High similarity to known drug'
                }
            
            return {
                **ttd_results,
                'warning': False,
                'message': 'No significant similarity found'
            }

        except Exception as e:
            logger.error(f"Similarity check failed for {smiles}: {str(e)}")
            return {'error': str(e)}

    def generate_fingerprint(self, smiles: str) -> Optional[DataStructs.cDataStructs.ExplicitBitVect]:
        """
        Generate 2048-bit Morgan fingerprint with radius 2
        Returns: RDKit fingerprint vector or None
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol or not mol.GetNumAtoms():
                return None
                
            return rdMolDescriptors.GetMorganFingerprintAsBitVect(
                mol, 
                radius=2,
                nBits=2048,
                useChirality=True
            )
        except Exception as e:
            logger.error(f"Fingerprint generation failed for {smiles}: {str(e)}")
            return None

    def query_ttd(self, fp: DataStructs.cDataStructs.ExplicitBitVect) -> Dict:
        """
        Query Therapeutic Target Database API with fingerprint
        Returns: Dict with similarity results or error
        """
        try:
            # Convert to on-bit indices for efficient transfer
            on_bits = fp.GetOnBits()
            
            response = requests.post(
                f"{self.ttd_endpoint}/similarity",
                json={
                    'on_bits': list(on_bits),
                    'n_bits': 2048,
                    'threshold': self.similarity_threshold
                },
                headers={'Authorization': f"Bearer {settings.TTD_API_KEY}"},
                timeout=15
            )
            response.raise_for_status()
            data = response.json()
            
            # Validate API response structure
            if not all(key in data for key in ['max_similarity', 'best_match']):
                raise ValueError("Invalid TTD response structure")
                
            best_match = data['best_match']
            if not best_match.get('drug_id'):
                raise ValueError("Missing drug ID in best match")

            return {
                'max_similarity': data['max_similarity'],
                'best_match_id': best_match['drug_id'],
                'best_match_name': best_match.get('name', 'N/A'),
                'mechanism_of_action': best_match.get('mechanism', 'Unknown')
            }
            
        except requests.JSONDecodeError as e:
            logger.error(f"TTD API JSON error: {str(e)}")
            return {'error': 'Invalid API response format'}
        except Exception as e:
            logger.error(f"TTD API request failed: {str(e)}")
            return {'error': f"TTD service error: {str(e)}"}

    def scrape_covdrug(self, drug_id: str) -> Dict:
        """
        Scrape CoviDrug website for repurposing information
        Returns: Dict with drug data or error
        """
        try:
            response = requests.get(
                f"{self.covidrug_base}/drug/{drug_id}",
                headers={'User-Agent': 'Mozilla/5.0'},
                timeout=15
            )
            response.raise_for_status()
            
            soup = BeautifulSoup(response.text, 'html.parser')
            
            return {
                'repurposing_status': self._parse_repurposing_status(soup),
                'clinical_trials': self._parse_clinical_trials(soup),
                'source_url': response.url
            }
            
        except requests.exceptions.RequestException as e:
            logger.error(f"CoviDrug request failed for {drug_id}: {str(e)}")
            return {'error': f"Network error: {str(e)}"}
        except Exception as e:
            logger.error(f"CoviDrug parsing failed for {drug_id}: {str(e)}")
            return {'error': 'Failed to parse drug information'}

    def _parse_repurposing_status(self, soup: BeautifulSoup) -> str:
        """Extract repurposing status from HTML"""
        status_div = soup.find('div', class_='repurposing-status')
        return status_div.text.strip() if status_div else 'Unknown'

    def _parse_clinical_trials(self, soup: BeautifulSoup) -> list:
        """Extract clinical trial information from HTML"""
        trials = []
        for card in soup.select('.trial-card'):
            try:
                trials.append({
                    'phase': card.select_one('.trial-phase').text,
                    'status': card.select_one('.trial-status').text,
                    'identifier': card.select_one('.trial-id').text
                })
            except AttributeError:
                continue
        return trials

    @staticmethod
    def calculate_tanimoto(fp1, fp2) -> float:
        """
        Calculate Tanimoto similarity between fingerprints
        Uses RDKit's optimized implementation
        """
        return DataStructs.TanimotoSimilarity(fp1, fp2)
