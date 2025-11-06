import streamlit as st

import modules.construct_design as construct_design

st.title("Design primers")

if "target_data" not in st.session_state:
    st.write(
        "No target data found, please go to the home page to retrieve data on a target protein."
    )

else:
    st.write("Designing primers for target: " + st.session_state.target_data.uniprot_id)
    st.markdown("### Primer table:")
    df = construct_design.generate_primer_dataframe(
        construct_dictionary=st.session_state.constructs,
        target_data=st.session_state.target_data,
    )
    st.session_state.primer_dataframe = df
    st.dataframe(data=df)
