use std::ptr::NonNull;

use super::bindings::tsk_id_t;
use super::bindings::tsk_site_table_add_row;
use super::bindings::tsk_site_table_clear;
use super::bindings::tsk_site_table_init;
use super::bindings::tsk_site_table_t;
use super::tskbox::TskBox;
use super::TskitError;

#[derive(Debug)]
pub struct SiteTable(TskBox<tsk_site_table_t>);

impl SiteTable {
    pub fn new(options: u32) -> Result<Self, TskitError> {
        let tsk =
            TskBox::new(|e: *mut tsk_site_table_t| unsafe { tsk_site_table_init(e, options) })?;
        Ok(Self(tsk))
    }

    pub unsafe fn new_borrowed(ptr: NonNull<tsk_site_table_t>) -> Self {
        let tsk = TskBox::new_init_from_ptr(ptr);
        Self(tsk)
    }

    pub fn as_ref(&self) -> &tsk_site_table_t {
        self.0.as_ref()
    }

    pub fn as_mut(&mut self) -> &mut tsk_site_table_t {
        self.0.as_mut()
    }

    pub fn clear(&mut self) -> i32 {
        unsafe { tsk_site_table_clear(self.as_mut()) }
    }

    pub fn add_row(
        &mut self,
        position: f64,
        ancestral_state: Option<&[u8]>,
    ) -> Result<tsk_id_t, TskitError> {
        self.add_row_with_metadata(position, ancestral_state, &[])
    }

    pub fn add_row_with_metadata(
        &mut self,
        position: f64,
        ancestral_state: Option<&[u8]>,
        metadata: &[u8],
    ) -> Result<tsk_id_t, TskitError> {
        unsafe {
            Ok(tsk_site_table_add_row(
                self.as_mut(),
                position,
                match ancestral_state {
                    Some(d) => d.as_ptr() as *const i8,
                    None => std::ptr::null(),
                },
                match ancestral_state {
                    Some(d) => d.len() as u64,
                    None => 0,
                },
                metadata.as_ptr().cast::<i8>(),
                metadata.len() as u64,
            ))
        }
    }
}

impl Default for SiteTable {
    fn default() -> Self {
        Self::new(0).unwrap()
    }
}
