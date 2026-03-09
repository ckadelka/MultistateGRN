#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# r"""
# This module provides functionality for retrieving, parsing, and loading
# biological Boolean network models from public online repositories.

# The :mod:`boolforge.bio_models` module allows users to programmatically access
# and import published Boolean and logical gene regulatory network models from
# GitHub repositories such as:

# - expert-curated (ckadelka): manually curated models from the
#   Design Principles of Gene Regulatory Networks repository.

# - pystablemotifs (jcrozum): models accompanying the PyStableMotifs library.

# - biodivine (sybila): models from the Sybila Biodivine Boolean Models repository.

# Functions are provided to:

# - Recursively list and download files from GitHub folders using the REST API.
# - Fetch raw text or byte content from remote sources.
# - Parse Boolean network models into BooleanNetwork objects.
# - Batch-download and convert all models from supported repositories.

# This module is intended to facilitate reproducible research by providing
# direct access to real-world Boolean GRN models for simulation, comparison,
# and benchmarking.

# Example
# -------
# >>> from boolforge import bio_models
# >>> result = bio_models.get_bio_models_from_repository('')
# >>> len(result['BooleanNetworks'])
# 122
# """

# import pickle
# import io

# try:
#     import requests
# except ImportError as e:
#     raise ImportError(
#         "The optional dependency 'requests' is required for bio_models. "
#         "Install it with `pip install requests`."
#     ) from e

# from .boolean_network import BooleanNetwork

# __all__ = [
#     'load_model',
#     'get_bio_models_from_repository',
# ]

# def _get_content_in_remote_folder(
#     url: str,
#     file_names: list,
#     file_download_urls: list
# ) -> None:
#     """
#     Recursively collect file names and raw download URLs from a GitHub folder.

#     Parameters
#     ----------
#     url : str
#         GitHub API URL pointing to a repository folder.
#     file_names : list
#         List that will be populated with discovered file names.
#     file_download_urls : list
#         List that will be populated with corresponding raw download URLs.

#     Returns
#     -------
#     None
#     """
#     import logging

#     folder = requests.get(url)
#     folder.raise_for_status()
#     folder_json = folder.json()

#     for item in folder_json:
#         if item['size'] > 0 and item['download_url'] is not None:
#             file_names.append(item['name'])
#             file_download_urls.append(item['download_url'])
#         else:
#             try:
#                 _get_content_in_remote_folder(
#                     item['url'], file_names, file_download_urls
#                 )
#             except Exception as e:
#                 logging.warning(
#                     "Failed to access subfolder at %s: %s",
#                     item.get('url', '<unknown>'),
#                     e,
#                 )

# def get_content_in_remote_folder(url: str) -> tuple:
#     """
#     Retrieve file names and raw download URLs from a GitHub repository folder.

#     Parameters
#     ----------
#     url : str
#         GitHub API URL pointing to a repository folder.

#     Returns
#     -------
#     tuple[list[str], list[str]]
#         A tuple ``(file_names, file_download_urls)``, where:

#         - ``file_names`` contains the names of discovered files.
#         - ``file_download_urls`` contains corresponding raw download URLs.
#     """
#     file_names = []
#     file_download_urls = []
#     _get_content_in_remote_folder(url, file_names, file_download_urls)
#     return file_names, file_download_urls


# def fetch_file(download_url: str) -> str:
#     """
#     Download raw text content from a remote file.

#     Parameters
#     ----------
#     download_url : str
#         Direct download URL to the file.

#     Returns
#     -------
#     str
#         File content as plain text.
#     """
#     r = requests.get(download_url)
#     r.raise_for_status()
#     return r.text


# def fetch_file_bytes(download_url: str) -> bytes:
#     """
#     Download raw binary content from a remote file.

#     Parameters
#     ----------
#     download_url : str
#         Direct download URL to the file.

#     Returns
#     -------
#     bytes
#         File content as raw bytes.
#     """
#     r = requests.get(download_url)
#     r.raise_for_status()
#     return r.content


# def load_model(
#     download_url: str,
#     max_degree: int = 24,
#     possible_separators: list = ['* =', '*=', '=', ','],
#     original_not: str = 'NOT',
#     original_and: str = 'AND',
#     original_or: str = 'OR',
#     ignore_first_line: bool = False,
#     simplify_functions: bool = False,
# ) -> BooleanNetwork:
#     """
#     Load and parse a Boolean network model from a remote text file.

#     Parameters
#     ----------
#     download_url : str
#         Direct download URL to the model file.
#     max_degree : int, optional
#         Maximum allowed in-degree for nodes (default: 24).
#     possible_separators : list[str], optional
#         Possible assignment separators used in the model file.
#     original_not : str, optional
#         Logical negation operator used in the model file.
#     original_and : str, optional
#         Logical AND operator used in the model file.
#     original_or : str, optional
#         Logical OR operator used in the model file.
#     ignore_first_line : bool, optional
#         If True, skip the first line of the file (default: False).
#     simplify_functions : bool, optional
#         If True, Boolean update functions are simplified after initialization.
#         Default is False.
        
#     Returns
#     -------
#     BooleanNetwork
#         Parsed Boolean network.

#     Raises
#     ------
#     ValueError
#         If the model cannot be parsed.
#     """
#     string = fetch_file(download_url)

#     if ignore_first_line:
#         string = string[string.index('\n') + 1:]

#     try:
#         bn = BooleanNetwork.from_string(
#             string,
#             possible_separators,
#             max_degree,
#             original_not,
#             original_and,
#             original_or,
#             simplify_functions=simplify_functions
#         )
#     except Exception as e:
#         raise ValueError(
#             f"Failed to parse Boolean network model from {download_url}"
#         ) from e

#     return bn

# def get_bio_models_from_repository(
#     repository: str = 'expert-curated (ckadelka)',
#     download_urls_pystablemotifs: list[str] | None = None,
#     max_degree: int = 24,
#     simplify_functions: bool = False,
# ) -> dict:
#     """
#     Load Boolean network models from selected online repositories.

#     This function downloads, parses, and constructs Boolean network models
#     from several curated online repositories. Models that cannot be parsed
#     are skipped and recorded separately.

#     Parameters
#     ----------
#     repository : str, optional
#         Identifier of the source repository. Supported values are:

#         - 'expert-curated (ckadelka)' (default)
#         - 'pystablemotifs (jcrozum)'
#         - 'biodivine (sybila)'

#     download_urls_pystablemotifs : list[str] or None, optional
#         Optional list of direct download URLs for PyStableMotifs models.
#         If provided, these URLs are used instead of querying the GitHub API (faster).
#         If None (default), model URLs are fetched dynamically from GitHub.
#     max_degree : int, optional
#         Maximum allowed in-degree for nodes (default: 24).
#     simplify_functions : bool, optional
#         If True, Boolean update functions are simplified after initialization.
#         Default is False.

#     Returns
#     -------
#     dict
#         Dictionary with the following keys:

#         - 'BooleanNetworks' : list[BooleanNetwork]
#             List of successfully parsed Boolean network models.

#         - 'SuccessfulDownloadURLs' : list[str]
#             URLs corresponding to models that were successfully loaded.

#         - 'FailedDownloadURLs' : list[str]
#             URLs corresponding to models that could not be parsed or loaded.
#     """
#     repositories = [
#         'expert-curated (ckadelka)',
#         'pystablemotifs (jcrozum)',
#         'biodivine (sybila)',
#     ]

#     bns = []
#     successful_download_urls = []
#     failed_download_urls = []

#     if repository == 'expert-curated (ckadelka)':
#         download_url_base = (
#             'https://raw.githubusercontent.com/ckadelka/'
#             'DesignPrinciplesGeneNetworks/main/'
#             'update_rules_122_models_Kadelka_SciAdv/'
#         )
#         download_url = download_url_base + 'all_txt_files.csv'
#         csv = fetch_file(download_url)

#         for line in csv.splitlines():
#             download_url = download_url_base + line
#             if '.txt' in download_url:
#                 try:
#                     if 'tabular' in download_url:
#                         F, I, var, constants = pickle.load(
#                             io.BytesIO(fetch_file_bytes(download_url))
#                         )
#                         for i in range(len(constants)):
#                             F.append([0, 1])
#                             I.append([len(var) + i])
#                         bn = BooleanNetwork(F, 
#                                             I, 
#                                             var + constants,
#                                             simplify_functions=simplify_functions)
#                     else:
#                         bn = load_model(
#                             download_url,
#                             original_and=" AND ",
#                             original_or=" OR ",
#                             original_not=" NOT ",
#                             simplify_functions=simplify_functions,
#                             max_degree = max_degree,
#                         )

#                     successful_download_urls.append(download_url)
#                     bns.append(bn)

#                 except Exception:
#                     failed_download_urls.append(download_url)

#     elif repository == 'pystablemotifs (jcrozum)':
#         if download_urls_pystablemotifs is None:
#             url = "https://api.github.com/repos/jcrozum/pystablemotifs/contents/models"
#             _, download_urls = get_content_in_remote_folder(url)
#         else:
#             download_urls = download_urls_pystablemotifs

#         for download_url in download_urls:
#             if '.txt' in download_url:
#                 try:
#                     bn = load_model(
#                         download_url,
#                         possible_separators=['*    =', '*   =', '*  =', '* =', '*='],
#                         original_and=[" and ", "&"],
#                         original_or=[" or ", "|"],
#                         original_not=[" not ", " !"],
#                         simplify_functions=simplify_functions,
#                         max_degree = max_degree,
#                     )
#                     successful_download_urls.append(download_url)
#                     bns.append(bn)
#                 except Exception:
#                     failed_download_urls.append(download_url)

#     elif repository == 'biodivine (sybila)':
#         download_url_base = (
#             'https://raw.githubusercontent.com/sybila/'
#             'biodivine-boolean-models/main/models/'
#         )
#         download_url = download_url_base + 'summary.csv'
#         csv = fetch_file(download_url)

#         for line in csv.splitlines():
#             try:
#                 ID, name, variables, inputs, regulations = line.split(', ')
#                 download_url = (
#                     download_url_base
#                     + '[id-%s]__[var-%s]__[in-%s]__[%s]/model.bnet'
#                     % (ID, variables, inputs, name)
#                 )
#                 bn = load_model(
#                     download_url,
#                     original_and=" & ",
#                     original_or=" | ",
#                     original_not="!",
#                     ignore_first_line=True,
#                     simplify_functions=simplify_functions,
#                     max_degree=max_degree,
#                 )
#                 successful_download_urls.append(download_url)
#                 bns.append(bn)
#             except Exception:
#                 failed_download_urls.append(download_url)

#     else:
#         raise ValueError(
#             "repository must be one of:\n - " + "\n - ".join(repositories)
#         )

#     return {
#         "BooleanNetworks": bns,
#         "SuccessfulDownloadURLs": successful_download_urls,
#         "FailedDownloadURLs": failed_download_urls,
#     }

