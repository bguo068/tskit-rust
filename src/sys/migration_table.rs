use std::ptr::NonNull;

use super::bindings::tsk_id_t;
use super::bindings::tsk_migration_table_add_row;
use super::bindings::tsk_migration_table_clear;
use super::bindings::tsk_migration_table_init;
use super::bindings::tsk_migration_table_t;
use super::tskbox::TskBox;
use super::Error;

#[derive(Debug)]
pub struct MigrationTable(TskBox<tsk_migration_table_t>);

impl MigrationTable {
    pub fn new(options: u32) -> Result<Self, Error> {
        let tsk = TskBox::new(|e: *mut tsk_migration_table_t| unsafe {
            tsk_migration_table_init(e, options)
        })?;
        Ok(Self(tsk))
    }

    pub unsafe fn new_borrowed(ptr: NonNull<tsk_migration_table_t>) -> Self {
        let tsk = TskBox::new_init_from_ptr(ptr);
        Self(tsk)
    }

    pub fn as_ref(&self) -> &tsk_migration_table_t {
        self.0.as_ref()
    }

    pub fn as_mut(&mut self) -> &mut tsk_migration_table_t {
        self.0.as_mut()
    }

    pub fn clear(&mut self) -> i32 {
        unsafe { tsk_migration_table_clear(self.as_mut()) }
    }

    pub fn add_row(
        &mut self,
        span: (f64, f64),
        node: tsk_id_t,
        source: tsk_id_t,
        dest: tsk_id_t,
        time: f64,
    ) -> Result<tsk_id_t, Error> {
        self.add_row_with_metadata(span, node, source, dest, time, &[])
    }

    pub fn add_row_with_metadata(
        &mut self,
        span: (f64, f64),
        node: tsk_id_t,
        source: tsk_id_t,
        dest: tsk_id_t,
        time: f64,
        metadata: &[u8],
    ) -> Result<tsk_id_t, Error> {
        unsafe {
            Ok(tsk_migration_table_add_row(
                self.as_mut(),
                span.0,
                span.1,
                node,
                source,
                dest,
                time,
                metadata.as_ptr().cast::<i8>(),
                metadata.len() as u64,
            ))
        }
    }
}

impl Default for MigrationTable {
    fn default() -> Self {
        Self::new(0).unwrap()
    }
}
