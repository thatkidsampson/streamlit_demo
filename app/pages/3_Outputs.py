import streamlit as st
from io import BytesIO
import pandas as pd
from streamlit_bokeh import streamlit_bokeh
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, FactorRange
from bokeh.transform import factor_cmap

import modules.construct_design as construct_design

st.title("Outputs")

if "target_data" not in st.session_state:
    st.write(
        "No target data found, please go to the home page to retrieve data on a target protein."
    )

else:
    st.markdown("### Primer order form:")
    primer_order_dataframe = construct_design.make_primer_plate(
        construct_df=st.session_state.construct_dataframe
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
        construct_df=st.session_state.construct_dataframe,
        primer_df=primer_order_dataframe,
    )
    st.dataframe(data=echo_input_dataframe)
    # Convert to a csv and show a download button
    csv = echo_input_dataframe.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download Echo input csv",
        csv,
        "echo_input.csv",
        "text/csv",
        key="download-echo-file",
    )

    # Plot the construct plate layout
    # assemble plot data
    plate_layout = pd.DataFrame(
        construct_design.generate_96_platemap(), columns=["Plate_well"]
    )
    plate_layout["row"] = plate_layout.Plate_well.str[:1]
    plate_layout["column"] = plate_layout.Plate_well.str[1:].astype(
        str
    )  # Keep as string throughout
    construct_data = st.session_state.construct_dataframe.copy()
    construct_data["Construct_name"] = construct_data.index
    plate_layout = pd.merge(plate_layout, construct_data, on="Plate_well", how="left")

    # Add a color column that's string-based for factor_cmap
    plate_layout["well_status"] = plate_layout["Construct_name"].apply(
        lambda x: "Filled" if pd.notna(x) else "Empty"
    )
    source = ColumnDataSource(plate_layout)

    # Create Bokeh figure
    st.markdown("### Construct plate layout:")
    construct_plate_layout = figure(
        height=450,
        width=670,
        title=None,
        x_axis_label="Column",
        y_axis_label="Row",
        toolbar_location=None,
        x_range=FactorRange(
            factors=list(map(str, sorted(plate_layout["column"].unique())))
        ),
        y_range=FactorRange(
            factors=list(reversed(sorted(plate_layout["row"].unique())))
        ),  # Reverse rows to match plate convention
    )

    # Plot with categorical coloring (filled vs empty wells)
    construct_plate_layout.scatter(
        x="column",
        y="row",
        source=source,
        size=30,
        fill_color=factor_cmap(
            "well_status", ["#e6e6e6", "#2ca02c"], ["Empty", "Filled"]
        ),
        line_color="black",
    )
    # Add hover tool to show well and construct info
    hover = HoverTool(
        tooltips=[("Plate well", "@Plate_well"), ("Construct name", "@Construct_name")]
    )
    construct_plate_layout.add_tools(hover)
    # Render Streamlit plot
    streamlit_bokeh(
        construct_plate_layout,
        use_container_width=False,
        theme="light_minimal",
        key="construct_plate_layout",
    )

    st.page_link("Home.py", label="Go back to construct design page", icon="💡")
    st.page_link("pages/2_Design_primers.py", label="Go to primer design page", icon="🧬")
