import streamlit as st
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
    primer_order_dataframe = construct_design.make_primer_plate(
        construct_df=st.session_state.primer_dataframe
    )
    st.dataframe(data=primer_order_dataframe)

    # Create a download button for the primer order form as an Excel file
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine="xlsxwriter")
    primer_order_dataframe.to_excel(writer, index=False)
    writer.close()
    output.seek(0)
    st.download_button(
        label="Download Primer Order Form",
        data=output.getvalue(),
        file_name="primer_order_form.xlsx",
        mime="application/vnd.ms-excel",
    )

    # Generate the Echo input file
    st.markdown("### Echo input file:")
    echo_input_dataframe = construct_design.make_echo_input_file(
        construct_df=st.session_state.primer_dataframe,
        primer_df=primer_order_dataframe)
    st.dataframe(data=echo_input_dataframe)
    # Convert to a csv and show a download button
    csv = echo_input_dataframe.to_csv(index=False).encode('utf-8')
    st.download_button(
        "Download Echo input csv",
        csv,
        "echo_input.csv",
        "text/csv",
        key='download-echo-file'
    )
