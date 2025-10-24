import streamlit as st
from streamlit_molstar import st_molstar_remote
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
from streamlit_bokeh import streamlit_bokeh

import modules.construct_design as construct_design

st.title("Construct design tool")

target_sequence = None

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
                target_data = construct_design.fetch_target_data(uniprot_id_input)
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

    with st.form(key="Register sequence boundaries"):
        submit_button = st.form_submit_button("Submit")
        if submit_button:
            if selection[0] not in st.session_state.N_term_boundaries:
                st.session_state.N_term_boundaries.append(selection[0])
                st.success(
                    f"Added N-terminal boundary: {sequence_list[selection[0]]}{selection[0] + 1}"
                )
            if selection[1] not in st.session_state.C_term_boundaries:
                st.session_state.C_term_boundaries.append(selection[1])
                st.success(
                    f"Added C-terminal boundary: {sequence_list[selection[1]]}{selection[1] + 1}"
                )
        if st.session_state.N_term_boundaries and st.session_state.C_term_boundaries:
            st.markdown(
                f"Current N-terminal boundaries: {[x + 1 for x in st.session_state.N_term_boundaries]}"
            )
            st.markdown(
                f"Current C-terminal boundaries: {[x + 1 for x in st.session_state.C_term_boundaries]}"
            )
        st.markdown(
            "Current number of constructs: "
            + str(
                len(st.session_state.N_term_boundaries)
                * len(st.session_state.C_term_boundaries)
            )
        )
        # assign construct names and assemble into a dictionary
        construct_number = 1
        sequence_length = st.session_state.target_data.sequence_length
        st.session_state.constructs = {
            st.session_state.target_data.uniprot_id: (1, sequence_length)
        }
        for Nterm in st.session_state.N_term_boundaries:
            for Cterm in st.session_state.C_term_boundaries:
                construct_name = f"{st.session_state.target_data.uniprot_id}_construct_{construct_number}"
                st.session_state.constructs[construct_name] = (Nterm + 1, Cterm + 1)
                construct_number += 1

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

    reset = st.button("Clear stored boundaries", on_click=clear_boundaries)
