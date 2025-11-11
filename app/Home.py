import streamlit as st
from streamlit_molstar import st_molstar_remote
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
from streamlit_bokeh import streamlit_bokeh

import modules.construct_design as construct_design

st.title("Construct design tool")
with st.expander("About this tool..."):
    st.write("""
        This tool is designed to help with the design of protein constructs for expression. \n
        Start by entering a UniProt ID to fetch the target protein sequence and structure prediction from AlphaFoldDB. \n
        Next, select N- and C-terminal boundaries for your constructs using the slider.\n
        You can add multiple boundaries, and the tool will generate all possible constructs based on your selections.\n
        Finally, proceed to the primer design page to generate primers for your constructs, and review the outputs on the outputs page. \n
        The tool will output a primer order file in the format required by the Merck primer order system and a picklist to dispense the required primers using an Echo liquid handler.
    """)

target_sequence = None

# number of wells availabke in the construct plate
PLATE_CAPACITY = 96

if "N_term_boundaries" not in st.session_state:
    st.session_state.N_term_boundaries = []
if "C_term_boundaries" not in st.session_state:
    st.session_state.C_term_boundaries = []


def clear_boundaries():
    st.session_state.N_term_boundaries = []
    st.session_state.C_term_boundaries = []


if "target_data" not in st.session_state:
    uniprot_id_input = st.text_input(
        "Enter a UniProt ID to fetch target data from AlphaFoldDB:", ""
    )
    if st.button("Fetch target structure prediction"):
        if uniprot_id_input:
            try:
                target_data = construct_design.fetch_target_data(
                    uniprot_id=uniprot_id_input
                )
                st.session_state.target_data = target_data
                st.success(f"Fetched data for UniProt ID: {uniprot_id_input}")
                st_molstar_remote(
                    st.session_state.target_data.alphafold_db_url, height=600
                )
            except RuntimeError as e:
                st.error(str(e))
else:
    st.write(
        "Designing constructs for target: " + st.session_state.target_data.uniprot_id
    )
    st_molstar_remote(st.session_state.target_data.alphafold_db_url, height=600)

if "target_data" in st.session_state:
    target_sequence = st.session_state.target_data.uniprot_sequence
    sequence_list = list(target_sequence)
    sequence_length = len(sequence_list)
    selection = st.select_slider(
        "Select N- and C-terminal construct boundaries",
        range(sequence_length),
        value=(0, sequence_length - 1),
        format_func=(lambda x: f"{sequence_list[x]}{x + 1}"),
    )

    # Display the target sequence with the selected sequence highlighted
    start_idx, end_idx = selection
    before_selection = target_sequence[:start_idx]
    selected_portion = target_sequence[start_idx : end_idx + 1]
    after_selection = target_sequence[end_idx + 1 :]

    highlighted_sequence = f"{before_selection}<mark style='background-color: #ffeb3b; padding: 2px; font-weight: bold;'>{selected_portion}</mark>{after_selection}"
    st.markdown(highlighted_sequence, unsafe_allow_html=True)
    col1, col2, col3 = st.columns(3)
    with col1:
        submit_button = st.button("Add construct boundaries", type="secondary")
    with col2:
        reset = st.button(
            "Clear stored boundaries", type="primary", on_click=clear_boundaries
        )
    with col3:
        st.page_link(
            "pages/2_Design_primers.py", label="Go to primer design page", icon="🧬"
        )
    if submit_button:
        if selection[0] not in st.session_state.N_term_boundaries:
            st.session_state.N_term_boundaries.append(selection[0])
            st.toast(
                f"Added N-terminal boundary: {sequence_list[selection[0]]}{selection[0] + 1}",
                icon="✅",
            )
        if selection[1] not in st.session_state.C_term_boundaries:
            st.session_state.C_term_boundaries.append(selection[1])
            st.toast(
                f"Added C-terminal boundary: {sequence_list[selection[1]]}{selection[1] + 1}",
                icon="✅",
            )
    st.markdown("#### Current construct info:")
    with st.container(key="Register sequence boundaries", border=True):
        # Display useful info about current boundaries and plate capacity
        col1, col2 = st.columns(2)
        number_of_constructs = len(st.session_state.N_term_boundaries) * len(
            st.session_state.C_term_boundaries
        )
        progress = number_of_constructs / PLATE_CAPACITY
        progress_text = f"Current number of constructs: {str(number_of_constructs)} / {PLATE_CAPACITY}"
        with col1:
            st.info(
                f"Current N-terminal boundaries: {[x + 1 for x in st.session_state.N_term_boundaries]}"
            )
            st.warning(
                f"Current C-terminal boundaries: {[x + 1 for x in st.session_state.C_term_boundaries]}"
            )
        with col2:
            if number_of_constructs > PLATE_CAPACITY:
                st.error(progress_text, icon="🚨")
                my_bar = st.progress(100, text="Plate over capacity.")
            else:
                st.success(progress_text, icon="✅")
                my_bar = st.progress(progress, text="Plate capacity used:")

        # assign construct names and assemble into a dictionary
        st.session_state.constructs = construct_design.generate_construct_dictionary(
            n_term_boundaries=st.session_state.N_term_boundaries,
            c_term_boundaries=st.session_state.C_term_boundaries,
            target_data=st.session_state.target_data,
        )

        # assemble plot data
        x_data = []
        y_data = []
        for construct, residue_range in st.session_state.constructs.items():
            for position in range(residue_range[0], residue_range[1]):
                y_data.append(construct)
                x_data.append(str(position))
        x_range = [str(x) for x in range(1, sequence_length)]
        y_range = list(st.session_state.constructs.keys())
        plot_height = 60 + (10 * len(y_range))
        source = ColumnDataSource(dict(x=x_data, y=y_data))

        # Create Bokeh figure
        construct_plot = figure(
            height=plot_height,
            title="Constructs vs sequence",
            x_axis_label="Residue number",
            y_axis_label="Construct",
            x_range=x_range,
            y_range=y_range,
            toolbar_location=None,
            tools="box_zoom, reset",
        )
        construct_plot.rect(x="x", y="y", width=0.8, height=0.4, source=source)
        hover = HoverTool(tooltips=[("Construct", "@y"), ("Residue", "@x")])
        construct_plot.add_tools(hover)
        construct_plot.xaxis.visible = False

        # Render Streamlit plot
        streamlit_bokeh(
            construct_plot,
            use_container_width=True,
            theme="streamlit",
            key="my_unique_key",
        )
