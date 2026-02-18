"""Streamlit dashboard entrypoint."""

import streamlit as st


def main() -> None:
    st.set_page_config(page_title="Immune Cell Dashboard", layout="wide")
    st.title("Immune Cell Clinical Trial Dashboard")
    st.caption("Dashboard scaffold is ready.")


if __name__ == "__main__":
    main()
