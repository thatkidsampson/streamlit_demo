import streamlit as st
from streamlit_molstar import st_molstar_remote

st.title("Construct design tool")

EXAMPLE_SEQUENCE = "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE"

st_molstar_remote("https://files.rcsb.org/view/1LOL.cif", height=600)

if 'N_term_boundaries' not in st.session_state:
    st.session_state.N_term_boundaries = []
if 'C_term_boundaries' not in st.session_state:
    st.session_state.C_term_boundaries = []

# Move slider and highlighting outside the form for real-time updates
sequence_list = list(EXAMPLE_SEQUENCE)
sequence_length = len(sequence_list)
selection = st.select_slider('Select N- and C-terminal construct boundaries',range(sequence_length), value=(0, sequence_length-1), format_func=(lambda x: f"{sequence_list[x]}{x+1}"))

# Create highlighted sequence display that updates in real-time
start_idx, end_idx = selection
before_selection = EXAMPLE_SEQUENCE[:start_idx]
selected_portion = EXAMPLE_SEQUENCE[start_idx:end_idx+1]
after_selection = EXAMPLE_SEQUENCE[end_idx+1:]

highlighted_sequence = f"{before_selection}<mark style='background-color: #ffeb3b; padding: 2px; font-weight: bold;'>{selected_portion}</mark>{after_selection}"
st.markdown(highlighted_sequence, unsafe_allow_html=True)

reset = st.button("Clear stored boundaries")
if reset:
    st.session_state.N_term_boundaries = []
    st.session_state.C_term_boundaries = []

with st.form(key="Register sequence boundaries"):
    submit_button = st.form_submit_button("Submit")
    if submit_button:
        if selection[0] not in st.session_state.N_term_boundaries:
            st.session_state.N_term_boundaries.append(selection[0])
            st.success(f"Added N-terminal boundary: {sequence_list[selection[0]]}{selection[0]+1}")
        if selection[1] not in st.session_state.C_term_boundaries:
            st.session_state.C_term_boundaries.append(selection[1])
            st.success(f"Added C-terminal boundary: {sequence_list[selection[1]]}{selection[1]+1}")
    if st.session_state.N_term_boundaries and st.session_state.C_term_boundaries:
        st.markdown(f"Current N-terminal boundaries: {[x+1 for x in st.session_state.N_term_boundaries]}")
        st.markdown(f"Current C-terminal boundaries: {[x+1 for x in st.session_state.C_term_boundaries]}")
    st.markdown("Current number of constructs: " + str(len(st.session_state.N_term_boundaries)*len(st.session_state.C_term_boundaries)))
