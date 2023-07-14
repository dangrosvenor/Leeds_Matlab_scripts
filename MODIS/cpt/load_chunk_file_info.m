function [am3_filepaths_chunk,gcm_year_cases,gcm_month_cases,gcm_model_cases] = load_chunk_file_info(gcm_saved_filename_file)


load(gcm_saved_filename_file,'am3_filepaths_chunk','gcm_year_cases','gcm_month_cases','gcm_model_cases');