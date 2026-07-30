#!/usr/bin/env python
"""
Lightweight file archive using ZIP for random access to thousands of signature files 
"""

import zipfile
import io
import shutil
import sys
import argparse
import logging
import logging
from pathlib import Path
from typing import List, Union, Optional, Iterator, Tuple
from contextlib import contextmanager

class FileArchive:
    """
    A lightweight archive manager for storing and accessing thousands of small files.
    Uses ZIP format for efficient random access.
    """
    
    def __init__(self, archive_path: Union[str, Path], mode: str = 'r'):
        """
        Initialize the archive.
        
        Args:
            archive_path: Path to the archive file
            mode: 'r' for read, 'w' for write (overwrites), 'a' for append
        """
        self.archive_path = Path(archive_path)
        self.mode = mode
        self._zipfile = None
        
    def __enter__(self):
        """Context manager entry."""
        self.open()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
        
    def open(self):
        """Open the archive."""
        compression = zipfile.ZIP_DEFLATED if self.mode in ['w', 'a'] else zipfile.ZIP_STORED
        self._zipfile = zipfile.ZipFile(self.archive_path, self.mode, compression=compression)
        
    def close(self):
        """Close the archive."""
        if self._zipfile:
            self._zipfile.close()
            self._zipfile = None
            
    def add_file(self, file_path: Union[str, Path], archive_name: Optional[str] = None):
        """
        Add a single file to the archive.
        
        Args:
            file_path: Path to the file to add
            archive_name: Name to use in archive (defaults to original filename)
        """
        if self.mode not in ['w', 'a']:
            raise ValueError("Archive must be opened in write or append mode")
            
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
            
        archive_name = archive_name or file_path.name
        self._zipfile.write(file_path, archive_name)
        
    def add_files(self, file_paths: List[Union[str, Path]], 
                  preserve_structure: bool = False,
                  base_path: Optional[Union[str, Path]] = None,
                  callback: Optional[callable] = None):
        """
        Add multiple files to the archive.
        
        Args:
            file_paths: List of file paths to add
            preserve_structure: If True, preserve directory structure
            base_path: Base path to strip when preserving structure
            callback: Optional callback function(current, total, filename)
        """
        if self.mode not in ['w', 'a']:
            raise ValueError("Archive must be opened in write or append mode")
            
        total = len(file_paths)
        base_path = Path(base_path) if base_path else None
        
        for idx, file_path in enumerate(file_paths, 1):
            file_path = Path(file_path)
            
            if preserve_structure and base_path:
                archive_name = str(file_path.relative_to(base_path))
            elif preserve_structure:
                archive_name = str(file_path)
            else:
                archive_name = file_path.name
                
            self.add_file(file_path, archive_name)
            
            if callback:
                callback(idx, total, str(file_path))
                
    def add_directory(self, directory: Union[str, Path], 
                     pattern: str = '*',
                     recursive: bool = True,
                     callback: Optional[callable] = None):
        """
        Add all files from a directory to the archive.
        
        Args:
            directory: Directory path
            pattern: File pattern (e.g., '*.txt', '*.json')
            recursive: If True, include subdirectories
            callback: Optional callback function(current, total, filename)
        """
        directory = Path(directory)
        
        if recursive:
            files = list(directory.rglob(pattern))
        else:
            files = list(directory.glob(pattern))
            
        files = [f for f in files if f.is_file()]
        self.add_files(files, preserve_structure=True, base_path=directory, callback=callback)
        
    def read(self, filename: str) -> bytes:
        """
        Read a file from the archive.
        
        Args:
            filename: Name of the file in the archive
            
        Returns:
            File contents as bytes
        """
        if self.mode not in ['r', 'a']:
            raise ValueError("Archive must be opened in read or append mode")
            
        try:
            return self._zipfile.read(filename)
        except KeyError:
            raise FileNotFoundError(f"File not found in archive: {filename}")
            
    def read_text(self, filename: str, encoding: str = 'utf-8') -> str:
        """
        Read a text file from the archive.
        
        Args:
            filename: Name of the file in the archive
            encoding: Text encoding (default: utf-8)
            
        Returns:
            File contents as string
        """
        return self.read(filename).decode(encoding)
        
    def concat(self,
               filenames: List[str], 
               separator: bytes = b'', 
               skip_missing: bool = False) -> Tuple[bytes, List[str], List[str]]:
        """
        Concatenate multiple files from the archive.
        
        Args:
            filenames: List of filenames to concatenate
            separator: Optional separator between files
            skip_missing: If True, skip files not found in archive; if False, raise error
            
        Returns:
            Tuple of (concatenated content as bytes, list of processed files, list of skipped files)
        """
        if self.mode not in ['r', 'a']:
            raise ValueError("Archive must be opened in read or append mode")

        result = io.BytesIO()
        processed_files = []
        skipped_files = []
        needs_separator = False
        preserve_gzip_stream = None

        for filename in filenames:
            try:
                with self._zipfile.open(filename, 'r') as member:
                    header = member.read(2)
                    is_gzip_member = header == b'\x1f\x8b'

                    # Preserve gzip members as a raw multi-member gzip stream, which
                    # keeps the output valid in the same way `cat file1.gz file2.gz` does.
                    if preserve_gzip_stream is None:
                        preserve_gzip_stream = is_gzip_member

                    if needs_separator and separator and not preserve_gzip_stream:
                        result.write(separator)

                    if header:
                        result.write(header)
                    shutil.copyfileobj(member, result, length=1024 * 1024)

                processed_files.append(filename)
                needs_separator = True

            except KeyError:
                if skip_missing:
                    skipped_files.append(filename)
                else:
                    raise FileNotFoundError(f"File not found in archive: {filename}")
                    
        return result.getvalue(), processed_files, skipped_files
        
    def concat_text(self, filenames: List[str], 
                   separator: str = '', 
                   encoding: str = 'utf-8',
                   skip_missing: bool = False) -> Tuple[str, List[str], List[str]]:
        """
        Concatenate multiple text files from the archive.
        
        Args:
            filenames: List of filenames to concatenate
            separator: Optional separator between files
            encoding: Text encoding (default: utf-8)
            skip_missing: If True, skip files not found in archive; if False, raise error
            
        Returns:
            Tuple of (concatenated content as string, list of processed files, list of skipped files)
        """
        content, processed, skipped = self.concat(
            filenames, 
            separator.encode(encoding),
            skip_missing=skip_missing
        )
        return content.decode(encoding), processed, skipped
        
    def extract(self, filename: str, output_path: Union[str, Path]):
        """
        Extract a file from the archive.
        
        Args:
            filename: Name of the file in the archive
            output_path: Where to save the extracted file
        """
        content = self.read(filename)
        Path(output_path).write_bytes(content)
        
    def extract_all(self, output_dir: Union[str, Path], callback: Optional[callable] = None):
        """
        Extract all files from the archive.
        
        Args:
            output_dir: Directory to extract files to
            callback: Optional callback function(current, total, filename)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        files = self.list_files()
        total = len(files)
        
        for idx, filename in enumerate(files, 1):
            self._zipfile.extract(filename, output_dir)
            
            if callback:
                callback(idx, total, filename)
                
    def list_files(self) -> List[str]:
        """
        List all files in the archive.
        
        Returns:
            List of filenames
        """
        return self._zipfile.namelist()
        
    def exists(self, filename: str) -> bool:
        """
        Check if a file exists in the archive.
        
        Args:
            filename: Name of the file to check
            
        Returns:
            True if file exists, False otherwise
        """
        return filename in self._zipfile.namelist()
        
    def get_info(self, filename: str) -> dict:
        """
        Get information about a file in the archive.
        
        Args:
            filename: Name of the file
            
        Returns:
            Dictionary with file info (size, compressed_size, etc.)
        """
        try:
            info = self._zipfile.getinfo(filename)
            return {
                'filename': info.filename,
                'size': info.file_size,
                'compressed_size': info.compress_size,
                'date_time': info.date_time,
                'compress_type': info.compress_type
            }
        except KeyError:
            raise FileNotFoundError(f"File not found in archive: {filename}")
            
    def get_size(self) -> int:
        """
        Get the total size of the archive file.
        
        Returns:
            Archive size in bytes
        """
        return self.archive_path.stat().st_size
        
    def __len__(self) -> int:
        """Return the number of files in the archive."""
        return len(self.list_files())
        
    def __contains__(self, filename: str) -> bool:
        """Check if a file exists in the archive."""
        return self.exists(filename)


# Convenience functions for quick operations

def create_archive(archive_path: Union[str, Path], 
                  file_paths: List[Union[str, Path]],
                  **kwargs) -> FileArchive:
    """
    Create a new archive with the given files.
    
    Args:
        archive_path: Path for the new archive
        file_paths: List of files to add
        **kwargs: Additional arguments for add_files()
        
    Returns:
        Closed FileArchive instance
    """
    with FileArchive(archive_path, 'w') as archive:
        archive.add_files(file_paths, **kwargs)
    return FileArchive(archive_path, 'r')


def quick_read(archive_path: Union[str, Path], filename: str) -> bytes:
    """
    Quick read of a single file from an archive.
    
    Args:
        archive_path: Path to the archive
        filename: Name of file to read
        
    Returns:
        File contents as bytes
    """
    with FileArchive(archive_path, 'r') as archive:
        return archive.read(filename)


def quick_concat(archive_path: Union[str, Path], 
                filenames: List[str],
                output_path: Optional[Union[str, Path]] = None,
                separator: bytes = b'',
                skip_missing: bool = True) -> Tuple[bytes, List[str], List[str]]:
    """
    Quick concatenation of files from an archive.
    
    Args:
        archive_path: Path to the archive
        filenames: List of files to concatenate
        output_path: Optional path to save the result
        separator: Optional separator between files
        skip_missing: If True, skip files not found in archive
        
    Returns:
        Tuple of (concatenated content as bytes, list of processed files, list of skipped files)
    """
    with FileArchive(archive_path, 'r') as archive:
        result, processed, skipped = archive.concat(filenames, separator=separator, skip_missing=skip_missing)
        
    if output_path:
        Path(output_path).write_bytes(result)
        
    return result, processed, skipped


def archive_genomes(genome_files: List[Union[str, Path]], 
                   output_archive: Union[str, Path],
                   preserve_structure: bool = False,
                   base_path: Optional[Union[str, Path]] = None,
                   verbose: bool = True) -> FileArchive:
    """
    Archive a list of genome files into a single ZIP file.
    
    Args:
        genome_files: List of genome file paths to archive
        output_archive: Path for the output archive file
        preserve_structure: If True, preserve directory structure in archive
        base_path: Base path to strip when preserving structure
        verbose: If True, print progress information
        
    Returns:
        FileArchive instance (closed)
        
    Example:
        >>> genome_files = ['genome1.fasta', 'genome2.fasta', 'genome3.fasta']
        >>> archive_genomes(genome_files, 'genomes.zip')
    """
    genome_files = [Path(f) for f in genome_files]
    
    # Validate all files exist
    missing_files = [f for f in genome_files if not f.exists()]
    if missing_files:
        raise FileNotFoundError(f"The following files were not found: {missing_files}")
    
    # Create progress callback if verbose
    def progress(current, total, filename):
        if verbose:
            logging.debug(f"[{current}/{total}] Adding: {filename}")
    
    callback = progress if verbose else None
    
    if verbose:
        logging.debug(f"Creating archive: {output_archive}")
        logging.debug(f"Number of genome files: {len(genome_files)}")
        total_size = sum(f.stat().st_size for f in genome_files)
        logging.debug(f"Total size: {total_size / (1024*1024):.2f} MB")
        logging.debug("-" * 50)
    
    # Create the archive
    with FileArchive(output_archive, 'w') as archive:
        archive.add_files(
            genome_files,
            preserve_structure=preserve_structure,
            base_path=base_path,
            callback=callback
        )
    
    if verbose:
        archive_size = Path(output_archive).stat().st_size
        compression_ratio = (1 - archive_size / total_size) * 100 if total_size > 0 else 0
        logging.debug("-" * 50)
        logging.debug(f"Archive created successfully!")
        logging.debug(f"Archive size: {archive_size / (1024*1024):.2f} MB")
        logging.debug(f"Compression ratio: {compression_ratio:.1f}%")
        logging.debug(f"Files in archive: {len(genome_files)}")
    
    return FileArchive(output_archive, 'r')


# Command-line interface
def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Archive genome files into a ZIP file with random access support.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Archive multiple genome files
  python sig_archive.py -o genomes.zip genome1.fasta genome2.fasta genome3.fasta
  
  # Archive all FASTA files from a directory
  python sig_archive.py -o genomes.zip -d /path/to/genomes -p "*.fasta"
  
  # Archive files from a list
  python sig_archive.py -o genomes.zip -l genome_list.txt
  
  # Archive with preserved directory structure
  python sig_archive.py -o genomes.zip -d /data/genomes --preserve
  
  # List files in an archive
  python sig_archive.py -i genomes.zip --list
  
  # Extract files from an archive
  python sig_archive.py -i genomes.zip --extract extracted/
  
  # Concatenate specific files (skip missing files by default)
  python sig_archive.py -i genomes.zip --concat genome1.fasta genome2.fasta -c output.fasta
  
  # Concatenate files from a list
  python sig_archive.py -i genomes.zip --concat-list files_to_concat.txt -c output.fasta
  
  # Concatenate with separator
  python sig_archive.py -i genomes.zip --concat-list files.txt -c output.fasta --separator "\\n---\\n"
  
  # Fail on missing files instead of skipping
  python sig_archive.py -i genomes.zip --concat-list files.txt -c output.fasta --strict
        """
    )
    
    # Input options
    input_group = parser.add_argument_group('Input Options')
    input_group.add_argument('files', nargs='*', help='Genome files to archive')
    input_group.add_argument('-d', '--directory', help='Directory containing genome files')
    input_group.add_argument('-p', '--pattern', default='*', 
                           help='File pattern for directory mode (default: *)')
    input_group.add_argument('-l', '--list', dest='file_list',
                           help='Text file containing list of genome files (one per line)')
    input_group.add_argument('-r', '--recursive', action='store_true',
                           help='Recursively search directory')
    
    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('-o', '--output', help='Output archive file path')
    output_group.add_argument('--preserve', action='store_true',
                            help='Preserve directory structure in archive')
    output_group.add_argument('-b', '--base-path', 
                            help='Base path to strip when preserving structure')
    
    # Archive operations
    operation_group = parser.add_argument_group('Archive Operations')
    operation_group.add_argument('-i', '--input-archive', 
                               help='Input archive file for reading operations')
    operation_group.add_argument('--list-files', action='store_true',
                               help='List files in the archive')
    operation_group.add_argument('--extract', metavar='DIR',
                               help='Extract all files to directory')
    operation_group.add_argument('--concat', nargs='+', metavar='FILE',
                               help='Concatenate specified files from archive')
    operation_group.add_argument('--concat-list', metavar='FILE',
                               help='Text file containing list of files to concatenate (one per line)')
    operation_group.add_argument('-c', '--concat-output',
                               help='Output file for concatenation')
    operation_group.add_argument('--separator', default='',
                               help='Separator to use between concatenated files (supports \\n, \\t)')
    operation_group.add_argument('--strict', action='store_true',
                               help='Fail on missing files instead of skipping them (concat mode)')
    
    # Other options
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Suppress progress output')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    return parser.parse_args()


def read_file_list(file_path: str, filename_only: bool = False) -> List[str]:
    """
    Read a list of filenames from a text file.
    
    Args:
        file_path: Path to the file containing the list
        filename_only: If True, return only filenames without paths (True by default for conacting files of the archive)
    Returns:
        List of filenames (stripped of whitespace)
    """
    with open(file_path, 'r') as f:
        if filename_only:
            return [line.strip().split('/')[-1] for line in f if line.strip() and not line.strip().startswith('#')]
        else:
            return [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]


def unescape_separator(separator: str) -> str:
    """
    Unescape common escape sequences in separator string.
    
    Args:
        separator: Separator string that may contain escape sequences
        
    Returns:
        Unescaped separator string
    """
    return separator.replace('\\n', '\n').replace('\\t', '\t').replace('\\r', '\r')


def main():
    """Main function for command-line interface."""
    args = parse_arguments()
    verbose = not args.quiet
    
    try:
        # Read operation mode
        if args.input_archive:
            if args.list_files:
                # List files in archive
                with FileArchive(args.input_archive, 'r') as archive:
                    files = archive.list_files()
                    logging.debug(f"Files in {args.input_archive}:")
                    logging.debug(f"Total: {len(files)} files")
                    logging.debug("-" * 50)
                    for filename in files:
                        info = archive.get_info(filename)
                        size_kb = info['size'] / 1024
                        logging.debug(f"{filename:50s} {size_kb:>10.2f} KB")
                    
            elif args.extract:
                # Extract files
                if verbose:
                    logging.debug(f"Extracting files from {args.input_archive} to {args.extract}")
                
                def progress(current, total, filename):
                    if verbose:
                        logging.debug(f"[{current}/{total}] Extracting: {filename}")
                
                with FileArchive(args.input_archive, 'r') as archive:
                    archive.extract_all(args.extract, callback=progress)
                
                if verbose:
                    logging.debug("Extraction completed!")
                    
            elif args.concat or args.concat_list:
                # Concatenate files
                if not args.concat_output:
                    logging.debug("Error: --concat-output (-c) is required when using --concat or --concat-list")
                    sys.exit(1)
                
                # Collect files to concatenate
                files_to_concat = []
                
                if args.concat:
                    files_to_concat.extend(args.concat)
                
                if args.concat_list:
                    list_files = read_file_list(args.concat_list)
                    files_to_concat.extend(list_files)
                
                if not files_to_concat:
                    logging.debug("Error: No files specified for concatenation")
                    sys.exit(1)
                
                # Remove duplicates while preserving order
                seen = set()
                unique_files = []
                for f in files_to_concat:
                    if f not in seen:
                        seen.add(f)
                        unique_files.append(f)
                
                if verbose:
                    logging.debug(f"Concatenating files from {args.input_archive}")
                    logging.debug(f"Requested files: {len(unique_files)}")
                    logging.debug("-" * 50)
                
                # Process separator
                separator = unescape_separator(args.separator).encode('utf-8')
                
                # Perform concatenation
                skip_missing = not args.strict
                result, processed, skipped = quick_concat(
                    args.input_archive, 
                    unique_files, 
                    args.concat_output,
                    separator=separator,
                    skip_missing=skip_missing
                )
                
                if verbose:
                    logging.debug(f"Successfully concatenated: {len(processed)} files")
                    if processed:
                        logging.debug("Processed files:")
                        for i, f in enumerate(processed, 1):
                            logging.debug(f"  {i}. {f}")
                    
                    if skipped:
                        logging.debug(f"\nSkipped (not found in archive): {len(skipped)} files")
                        logging.debug("Skipped files:")
                        for i, f in enumerate(skipped, 1):
                            logging.debug(f"  {i}. {f}")
                    
                    output_size = len(result) / 1024
                    logging.debug("-" * 50)
                    logging.debug(f"Output size: {output_size:.2f} KB")
                    logging.debug(f"Saved to: {args.concat_output}")
                
                # Exit with error code if some files were skipped and user might want to know
                if skipped and not verbose:
                    logging.warning(f"Warning: {len(skipped)} file(s) were not found in the archive and were skipped")
                    
            else:
                logging.error("Error: Specify an operation (--list-files, --extract, --concat, or --concat-list)")
                sys.exit(1)
                
        # Write operation mode
        else:
            if not args.output:
                logging.error("Error: --output is required when creating an archive")
                sys.exit(1)
            
            # Collect genome files
            genome_files = []
            
            if args.files:
                genome_files.extend(args.files)
            
            if args.directory:
                directory = Path(args.directory)
                if args.recursive:
                    found_files = list(directory.rglob(args.pattern))
                else:
                    found_files = list(directory.glob(args.pattern))
                genome_files.extend([str(f) for f in found_files if f.is_file()])
            
            if args.file_list:
                list_files = read_file_list(args.file_list)
                genome_files.extend(list_files)
            
            if not genome_files:
                logging.error("Error: No genome files specified")
                logging.error("Use --help for usage information")
                sys.exit(1)
            
            # Remove duplicates while preserving order
            seen = set()
            unique_files = []
            for f in genome_files:
                if f not in seen:
                    seen.add(f)
                    unique_files.append(f)
            
            # Create archive
            archive_genomes(
                unique_files,
                args.output,
                preserve_structure=args.preserve,
                base_path=args.base_path,
                verbose=verbose
            )
            
    except Exception as e:
        logging.error(f"Error: {e}")
        import traceback
        if verbose:
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
