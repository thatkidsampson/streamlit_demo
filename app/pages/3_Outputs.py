import streamlit as st
import xlsxwriter
from io import BytesIO
import pandas as pd

import modules.construct_design as construct_design

st.title("Outputs")

if "target_data" not in st.session_state:
    st.write(
        "No target data found, please go to the home page to retrieve data on a target protein."
    )

else:
    st.markdown("### Primer order form:")
    primer_order_dataframe = construct_design.make_primer_plate(input_df=st.session_state.primer_dataframe)
    st.dataframe(data=primer_order_dataframe)

    # Create a download button for the primer order form as an Excel file
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    primer_order_dataframe.to_excel(writer, index=False)
    writer.close()
    output.seek(0)
    st.download_button(
        label="Download Primer Order Form",
        data=output.getvalue(),
        file_name="primer_order_form.xlsx",
        mime="application/vnd.ms-excel"
    )