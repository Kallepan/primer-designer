use rusqlite::{Connection, Result};
use std::error::Error;

#[allow(dead_code)]
pub fn load_sqlite_from_file(input_file_path: &str) -> Result<Connection, Box<dyn Error>> {
    let conn = Connection::open(input_file_path)?;
    Ok(conn)
}