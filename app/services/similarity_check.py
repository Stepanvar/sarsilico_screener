import logging
from typing import Dict
import requests
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from bs4 import BeautifulSoup, Tag
from django.conf import settings

logger = logging.getLogger(__name__)

def check_similarity(compound: str, is_smile: bool) -> Dict:
    """
    Check similarity for a compound.

    Parameters:
      compound (str): Either a SMILES string or a sequence.
      is_smile (bool): True if the compound is a SMILES string, False if it is a sequence.

    Returns a dictionary containing:
      - max_similarity: The Tanimoto similarity score (as float).
      - best_match_id: The drug ID from TTD.
      - best_match_name: The drug name.
      - mechanism_of_action: Currently 'Unknown'.
      - If the similarity is below the threshold (defined in settings),
        the dictionary also gets a 'medicine_name' key via lookup_medicine_name.
      - In case of errors, an 'error' key is returned.
    """
    try:
        # Validate input
        if is_smile:
            mol = Chem.MolFromSmiles(compound)
            if not mol or mol.GetNumAtoms() == 0:
                logger.error("Invalid SMILES structure provided.")
                return {'error': 'Invalid SMILES structure'}
            # Generate fingerprint for validation (if needed)
            _ = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, useChirality=True)
        else:
            if not compound.strip():
                logger.error("Empty sequence provided.")
                return {'error': 'Empty sequence provided'}

        # Query the TTD similarity search page
        similarity_results = _query_ttd_api(compound)
        if 'error' in similarity_results:
            return similarity_results

        # Convert the similarity value to float (defaulting to 0.0 if conversion fails)
        try:
            similarity_value = float(similarity_results.get('max_similarity', 0))
        except (ValueError, TypeError):
            similarity_value = 0.0

        similarity_results['max_similarity'] = similarity_value

        # If similarity is larger than the threshold, call lookup_medicine_name
        if similarity_value > settings.SIMILARITY_THRESHOLD:
            medicine_name = similarity_results.get('best_match_name')
            similarity_results['medicine_name'] = medicine_name

        # Set a warning flag if similarity meets or exceeds the threshold
        similarity_results['warning'] = similarity_value >= settings.SIMILARITY_THRESHOLD

        return similarity_results

    except Exception as e:
        logger.error(f"Similarity check failed: {str(e)}")
        return {'error': 'Processing error'}

def _query_ttd_api(compound: str) -> Dict:
    """
    Internal function to query the TTD similarity search page.

    Performs a GET request to obtain hidden form fields, then posts the compound.
    Parses the returned HTML table to extract the first result row data.
    Expects the result table to have a header and then data rows with class "panel4show".

    Returns a dictionary with keys:
      - max_similarity: Float value from the third cell of the first data row.
      - best_match_id: Drug ID (from the first cell, using an <a> if available).
      - best_match_name: Drug name (from the second cell).
      - mechanism_of_action: 'Unknown'
    
    In case of error, returns a dictionary with an 'error' key.
    """
    try:
        url = "https://db.idrblab.net/ttd/ttd-search/drug-similarity"
        headers = {'User-Agent': 'Mozilla/5.0'}

        # GET the page to capture hidden form fields.
        get_response = requests.get(url, headers=headers, timeout=10)
        get_response.raise_for_status()
        soup = BeautifulSoup(get_response.text, 'html.parser')

        form_build_id_input = soup.find("input", {"name": "form_build_id"})
        form_id_input = soup.find("input", {"name": "form_id"})
        if not (isinstance(form_build_id_input, Tag) and isinstance(form_id_input, Tag)):
            logger.error("Required form fields not found on TTD page.")
            return {'error': 'Required form fields not found'}

        form_build_id = form_build_id_input.get("value", "")
        form_id = form_id_input.get("value", "")

        # POST the compound.
        payload = {
            'form_build_id': form_build_id,
            'form_id': form_id,
            # TTD expects the field "smiles" even for sequences.
            'smiles': compound,
            'op': 'Submit'
        }
        post_response = requests.post(url, data=payload, headers=headers, timeout=50)
        post_response.raise_for_status()

        result_soup = BeautifulSoup(post_response.text, 'html.parser')
        result_table = result_soup.find("table", class_="ttd-table-intro")
        if not (result_table and isinstance(result_table, Tag)):
            logger.error("No similarity results found in TTD response.")
            return {'error': 'No similarity results found'}

        # Use select() on the table (we now ensure result_table is a Tag)
        data_rows = result_table.select("tr.panel4show")
        if not data_rows:
            logger.error("No similarity result rows found in TTD results.")
            return {'error': 'No similarity result rows found'}

        first_result = data_rows[0]
        # Ensure first_result is a Tag before calling find_all()
        if not isinstance(first_result, Tag):
            logger.error("First result row is not a valid Tag.")
            return {'error': 'Unexpected result format'}

        cells = first_result.find_all("td")
        if len(cells) < 3:
            logger.error("Incomplete similarity result data.")
            return {'error': 'Incomplete similarity result data'}

        # Extract best match ID from the first cell
        first_cell = cells[0]
        if isinstance(first_cell, Tag):
            a_tag = first_cell.find("a")
            if a_tag and isinstance(a_tag, Tag):
                best_match_id = a_tag.get_text(strip=True)
            else:
                best_match_id = first_cell.get_text(strip=True)
        else:
            best_match_id = str(first_cell).strip()

        # Extract best match name from the second cell.
        second_cell = cells[1]
        best_match_name = second_cell.get_text(strip=True) if isinstance(second_cell, Tag) else str(second_cell).strip()

        # Extract Tanimoto similarity from the third cell.
        third_cell = cells[2]
        try:
            cell_text = third_cell.get_text(strip=True) if isinstance(third_cell, Tag) else str(third_cell).strip()
            tanimoto_value = float(cell_text)
        except (ValueError, TypeError):
            tanimoto_value = 0.0

        return {
            'max_similarity': tanimoto_value,
            'best_match_id': best_match_id,
            'best_match_name': best_match_name,
            'mechanism_of_action': 'Unknown'
        }

    except requests.RequestException as e:
        logger.error(f"TTD API error: {str(e)}")
        return {'error': 'Database connection failed'}
    except Exception as e:
        logger.error(f"Error processing TTD response: {str(e)}")
        return {'error': 'Data format error'}

def lookup_medicine_name(similarity_results: Dict) -> str:
    """
    Scrape the drug name from the CoviDrug website using meta tags.

    Uses the best_match_id from the similarity results to construct the URL.
    """
    drug_id = similarity_results.get('best_match_name')
    if not drug_id:
        return ''
    try:
        url = f"{settings.COVIDRUG_SCRAPE_URL}/drug/{drug_id}"
        response = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'}, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        meta_tags = soup.find_all('meta', {'property': 'og:title'})
        if not meta_tags:
            return ""
        meta_tag = meta_tags[0]
        if isinstance(meta_tag, Tag):
            # Ensure the return type is a string. If the attribute value is a list, join it.
            content = meta_tag.get('content', '')
            if isinstance(content, list):
                return " ".join(content)
            return content or ""
        else:
            return str(meta_tag)
    except Exception as e:
        logger.warning(f"Name lookup failed for {drug_id}: {str(e)}")
        return ''
