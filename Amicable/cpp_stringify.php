<?php
if ($argc < 2)
{
	exit();
}

$input_name = $argv[1];
$output_name = "$input_name.inl";
if (file_exists($output_name) && (filemtime($input_name) <= filemtime($output_name)))
{
	exit();
}

$data = str_replace("\r", '', file_get_contents($input_name));
$checksum = crc32($data);
$data = str_replace("\\", "\\\\", $data);
$data = str_replace("\"", "\\\"", $data);
$data = str_replace("\n", "\\n\"\\\n\"", $data);
file_put_contents($output_name, "#pragma once\n\nstatic const char* ".str_replace('.', '_', $input_name)." = \"$data\";\n\nstatic const unsigned int ".str_replace('.', '_', $input_name).'_crc32 = 0x'.dechex($checksum).";\n");
?>
