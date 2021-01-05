using System;
using System.IO;
using Microsoft.VisualBasic.FileIO;
using System.Collections.Generic;
using System.Linq;
using System.Xml.Linq;
using System.Text;
using System.Globalization;
using System.Threading.Tasks;
using System.Reflection;

namespace sgRNA
{
	class Program
	{
		/// <summary>
		/// SnapGeneファイルのセグメントを読み込むためのクラス
		/// </summary>
		public class SGSegment
		{
			public Byte SGIdentifier { get; set; }
			public Int32 SGSize { get; set; }
			public Byte[] SGContent { get; set; }
			public Byte SGSeqType { get; set; }
			public Byte[] SGSeq { get; set; }

			public SGSegment() { }
			public SGSegment(byte identifier, int size, byte[] content)
			{
				SGIdentifier = identifier;
				SGSize = size;
				SGContent = content;
			}
			public SGSegment(byte identifier, int size, byte[] content, byte type, byte[] seq)
			{
				SGIdentifier = identifier;
				SGSize = size;
				SGContent = content;
				SGSeqType = type;
				SGSeq = seq;
			}
		}

		/// <summary>
		/// SnapGeneファイルからセグメント単位で読み込みます
		/// </summary>
		/// <param name="FileName">SnapGeneファイルパス</param>
		/// <returns>全セグメントを格納したリスト</returns>
		public static List<SGSegment> ReadSegments(string FileName)
		{
			// 全セグメントを格納するリスト
			List<SGSegment> output = new List<SGSegment>();

			using (FileStream fs = new FileStream(FileName, FileMode.Open, FileAccess.Read))
			using (BinaryReader br = new BinaryReader(fs))
			{
				int fileSize = (int)fs.Length;
				while (br.BaseStream.Position != fileSize)
				{
					byte id = br.ReadByte();
					byte[] sizeByteArray = br.ReadBytes(4);
					Array.Reverse(sizeByteArray);
					int size = BitConverter.ToInt32(sizeByteArray, 0);

					SGSegment seg = new SGSegment();

					if (id != 0)
					{
						byte[] content = br.ReadBytes(size);
						seg = new SGSegment(id, size, content);
					}
					else
					{
						byte[] content = br.ReadBytes(size);
						byte type = content[0];
						byte[] seq = new byte[content.Length - 1];
						Array.Copy(content, 1, seq, 0, content.Length - 1);
						seg = new SGSegment(id, size, content, type, seq);
					}

					output.Add(seg);
				}
			}

			return output;
		}

		/// <summary>
		/// SnapGeneファイルを作成します
		/// </summary>
		/// <param name="SGsegs">SnapGeneセグメントを格納したリスト</param>
		/// <param name="FileName">保存先のファイルパス</param>
		/// <returns>ファイルが作成されたらtrueを返します</returns>
		public static bool SaveSGFiles(List<SGSegment> SGsegs, string FileName)
		{
			if (SGsegs == null || FileName == "") return false;

			using (FileStream fs = new FileStream(FileName, FileMode.Create))
			using (BinaryWriter bw = new BinaryWriter(fs))
			{
				foreach (SGSegment seg in SGsegs)
				{
					bw.Write(seg.SGIdentifier);
					byte[] sizeByteArray = BitConverter.GetBytes(seg.SGSize); // LittleEndian
					Array.Reverse(sizeByteArray); // BigEndian
					bw.Write(sizeByteArray);
					bw.Write(seg.SGContent);
				}
			}
			return true;
		}

		/// <summary>
		/// DNA配列フラグを文字列に変換して返す
		/// </summary>
		/// <param name="type">DNA配列フラグを格納したバイト型変数</param>
		/// <returns>変換後の文字列</returns>
		public static string TypeString(byte type)
		{
			string output = "";
			if ((type & 1) != 0)
			{
				output += "Circular, ";
			}
			else
			{
				output += "Linear, ";
			}
			if ((type & 2) != 0)
			{
				output += "Double-stranded, ";
			}
			else
			{
				output += "Single-stranded, ";
			}
			if ((type & 4) != 0)
			{
				output += "dam-methylated, ";
			}
			else
			{
				output += "dam-non-methylated, ";
			}
			if ((type & 8) != 0)
			{
				output += "dcm-methylated, ";
			}
			else
			{
				output += "dcm-non-methylated, ";
			}
			if ((type & 16) != 0)
			{
				output += "ecoKI-methylated, ";
			}
			else
			{
				output += "ecoKI-non-methylated, ";
			}

			return output;
		}

		public static List<SGSegment> PrepareSGSegments(List<SGSegment> vec, string enz, string name, List<byte> pattern, string seq_s)
        {
			// ベクターの配列を取得
			List<byte> seq_original = new List<byte>();
			foreach (SGSegment seg in vec)
            {
				if (seg.SGIdentifier == 0)
                {
					seq_original = seg.SGContent.ToList<byte>();

					// センス鎖から5'オーバーハングを除去した配列を取得
					string seq_s_nooverhangs = seq_s.Remove(0, 4);
					byte[] seqBytes = System.Text.Encoding.UTF8.GetBytes(seq_s_nooverhangs);

					// ベクター配列内のNNN…NNNをオリゴ配列で置換したものを取得
					byte[] seq_replaced = Replace(seq_original.ToArray(), pattern.ToArray(), seqBytes);

					seg.SGContent = seq_replaced;
				}
            }

			return vec;
        }

		private static byte[] Replace(byte[] input, byte[] pattern, byte[] replacement)
        {
			if (pattern.Length == 0)
            {
				return input;
            }

			List<byte> result = new List<byte>();

			int i;
			for (i = 0; i <= input.Length - pattern.Length; i++)
			{
				bool foundMatch = true;
				for (int j = 0; j < pattern.Length; j++)
				{
					if (input[i + j] != pattern[j])
					{
						foundMatch = false;
						break;
					}
				}

				if (foundMatch)
				{
					result.AddRange(replacement);
					i += pattern.Length - 1;
				}
				else
				{
					result.Add(input[i]);
				}
			}
			
			for (; i < input.Length; i++)
            {
				result.Add(input[i]);
            }

			return result.ToArray();
        }

		private static void Main(string[] args)
		{
			Console.WriteLine("======================================================");
			Console.WriteLine("  CRISPR/Cas9用ツール");
			Console.WriteLine("    Copyright (c) 2018-2021 Genta Ito");
			Console.WriteLine("    Version 3.1");
			Console.WriteLine("======================================================");
			Console.WriteLine("  入力ファイル書式：");
			Console.WriteLine("    「ベクター名称+改行」");
			Console.WriteLine("    「標的配列名称1+タブ+標的配列1（5' -> 3'）+改行」");
			Console.WriteLine("    「標的配列名称2+タブ+標的配列2（5' -> 3'）+改行」");
			Console.WriteLine("    「ベクター名称+改行」");
			Console.WriteLine("    「…（以下同じ）」");
			Console.WriteLine("    各項目はタブで区切る。singleベクターの場合には名称2、配列2は省略可");
			Console.WriteLine("    対応singleベクター:pX335, pBabe Puro");
			Console.WriteLine("    対応dualベクター:pX459dual D10A");
			Console.WriteLine("======================================================");
			Console.WriteLine("  出力ファイル：");
			Console.WriteLine("  ・発注用オリゴDNAリスト（入力ファイルに追記）");
			Console.WriteLine("  ・それぞれのオリゴDNA配列SnapGeneファイル");
			Console.WriteLine("  ・アニーリングしたオリゴDNA配列SnapGeneファイル");
			Console.WriteLine("");
			Console.WriteLine("  入力ファイルがあるフォルダにファイルが作成されます");
			Console.WriteLine("======================================================");
			Console.WriteLine("");

			// ドラッグアンドドロップされたファイルのファイルパスを取得
			// 先頭に格納される実行ファイル名を除く
			string[] filePath = Environment.GetCommandLineArgs();
			int startIndex = 0;
			int numberOfFiles = 0;
			for (var i = 0; i < filePath.Length; i++)
			{
				int len = filePath[i].Length;
				if (filePath[i].Substring(len - 3, 3) != "exe" && filePath[i].Substring(len - 3, 3) != "dll")
				{
					startIndex = i;
					break;
				}
			}
			numberOfFiles = filePath.Length - startIndex;

			// エラー処理
			if (numberOfFiles != 1)
			{
				Console.WriteLine("  エラー：ドロップするファイルは1個だけにしてください。");
				Console.WriteLine("  終了するには何かキーを押してください。");
				Console.ReadKey();
				Environment.Exit(0);
			}

			// 入力ファイルのディレクトリを取得
			string dirName = Path.GetDirectoryName(filePath[startIndex]);
			// 実行ファイルのディレクトリを取得
			string exeDirName = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);

			// クエリディクショナリ
			var listNameSeq = new Dictionary<string, string>();
			// 結果ディクショナリ
			var oligos = new Dictionary<string, string>();
			var annealed_oligos = new Dictionary<string, string>();

			using (TextFieldParser parser = new TextFieldParser(filePath[startIndex], Encoding.ASCII))
			{
				try
				{
					parser.TextFieldType = FieldType.Delimited;
					parser.SetDelimiters("\t");

					while (parser.EndOfData == false)
					{
						string[] vector = parser.ReadFields();

						// RESULTSヘッダーがある場合はエラー
						if (vector[0] == "RESULTS")
                        {
							Console.WriteLine("  {0}は処理済みのファイルのようです。", filePath[startIndex]);
							Console.WriteLine("  終了するには何かキーを押してください。");
							Console.ReadKey();
							Environment.Exit(0);
						}
						else if (vector[0] == "pX335")
                        {
							string[] data = parser.ReadFields();
							listNameSeq.Add(data[0], data[1]);

							Dictionary<string, string> oligo = MakeOligoSeq(data[0], data[1]);
							oligos = oligos.Concat(oligo).ToDictionary(x => x.Key, x => x.Value);
							annealed_oligos = annealed_oligos.Concat(AnnealOligos(data[0], oligo[data[0] + "_s"])).ToDictionary(x => x.Key, x => x.Value);

							List<SGSegment> pX335 = new List<SGSegment>();

							List<byte> pattern = new List<byte>();

							if (oligo[data[0] + "_s"].Length == 24)
							{
								pX335 = ReadSegments(exeDirName + "\\pX335_N20.dna");
								pattern.AddRange(System.Text.Encoding.UTF8.GetBytes("NNNNNNNNNNNNNNNNNNNN").ToList<byte>()); // N20
							}
							else if (oligo[data[0] + "_s"].Length == 25)
                            {
								pX335 = ReadSegments(exeDirName + "\\pX335_N21.dna");
								pattern.AddRange(System.Text.Encoding.UTF8.GetBytes("NNNNNNNNNNNNNNNNNNNNN").ToList<byte>()); // N21
							}
							else { }

							List<SGSegment> SGsegs = PrepareSGSegments(pX335, "BbsI", data[0], pattern, oligo[data[0] + "_s"]);
							SaveSGFiles(SGsegs, dirName + "\\" + data[0] + "_pX335.dna");
						}
						else if (vector[0] == "pBabe Puro")
                        {
							string[] data = parser.ReadFields();
							listNameSeq.Add(data[0], data[1]);

							Dictionary<string, string> oligo = MakeOligoSeq(data[0], data[1]);
							oligos = oligos.Concat(oligo).ToDictionary(x => x.Key, x => x.Value);
							annealed_oligos = annealed_oligos.Concat(AnnealOligos(data[0], oligo[data[0] + "_s"])).ToDictionary(x => x.Key, x => x.Value);
						}
						else if (vector[0] == "pX459dual D10A")
                        {
							string[] data1 = parser.ReadFields();
							listNameSeq.Add(data1[0], data1[1]);

							Dictionary<string, string> oligo1 = MakeOligoSeq(data1[0], data1[1]);
							oligos = oligos.Concat(oligo1).ToDictionary(x => x.Key, x => x.Value);
							annealed_oligos = annealed_oligos.Concat(AnnealOligos(data1[0], oligo1[data1[0] + "_s"])).ToDictionary(x => x.Key, x => x.Value);

							string[] data2 = parser.ReadFields();
							listNameSeq.Add(data2[0], data2[1]);

							Dictionary<string, string> oligo2 = MakeOligoSeq(data2[0], data2[1]);
							oligos = oligos.Concat(oligo2).ToDictionary(x => x.Key, x => x.Value);
							annealed_oligos = annealed_oligos.Concat(AnnealOligos(data2[0], oligo2[data2[0] + "_s"])).ToDictionary(x => x.Key, x => x.Value);
						}
						else
                        {
							Console.WriteLine("  ベクター名称が不正です。");
							Console.WriteLine("  終了するには何かキーを押してください。");
							Console.ReadKey();
							Environment.Exit(0);
						}
					}
				}
				catch (Exception ex)
				{
					Console.WriteLine("  " + ex.Message);
					Console.WriteLine("  終了するには何かキーを押してください。");
					Console.ReadKey();
					Environment.Exit(0);
				}
			}

			// 結果ファイル作成
			var sw = new StreamWriter(filePath[startIndex], true);
			try
			{
				CultureInfo ci = CultureInfo.CurrentCulture;
				TextInfo ti = ci.TextInfo;
				DateTime dt = DateTime.Now;

				sw.WriteLine("");
				sw.WriteLine("");
				sw.WriteLine("RESULTS");
				sw.WriteLine(dt.ToString("yyyy/MM/dd HH:mm:ss"));
				sw.WriteLine("======== ユーロフィン発注用 ========");
				foreach (KeyValuePair<string, string> kv in oligos)
				{
					// 塩基配列はすべて大文字に変換
					sw.WriteLine("{0}" + "\t" + "{1}", kv.Key, ti.ToUpper(kv.Value));
				}
				sw.WriteLine("======== 以下余白 ========");

				Console.WriteLine("  入力ファイルに追記しました。");
			}
			catch (Exception ex)
			{
				Console.WriteLine("  " + ex.Message);
				Console.WriteLine("  終了するには何かキーを押してください。");
				Console.ReadKey();
				Environment.Exit(0);
			}
			finally
			{
				sw.Close();
			}

			// それぞれのオリゴDNA配列SnapGeneファイル作成
			foreach (KeyValuePair<string, string> kv in oligos)
			{
				try
				{
					string path = dirName + "\\" + kv.Key + ".dna";
					MakeSGFile(path, kv.Value, false);
					Console.WriteLine("  SnapGeneファイルを作成しました：{0}", path);
				}
				catch (Exception ex)
                {
					Console.WriteLine("  " + ex.Message);
					Console.WriteLine("  終了するには何かキーを押してください。");
					Console.ReadKey();
					Environment.Exit(0);
				}
			}

			// アニーリング済みオリゴDNA配列SnapGeneファイル作成
			foreach (KeyValuePair<string, string> kv in annealed_oligos)
            {
				try
                {
					string path = dirName + "\\" + kv.Key + ".dna";
					MakeSGFile(path, kv.Value, true);
					Console.WriteLine("  アニーリング済みオリゴDNA配列のSnapGeneファイルを作成しました：{0}", path);
				}
				catch (Exception ex)
				{
					Console.WriteLine("  " + ex.Message);
					Console.WriteLine("  終了するには何かキーを押してください。");
					Console.ReadKey();
					Environment.Exit(0);
				}
			}

			Console.WriteLine("  全ての処理が完了しました。");
			Console.WriteLine("  終了するには何かキーを押してください。");
			Console.ReadKey();
			Environment.Exit(0);
		}

		/// <summary>
		/// 標的配列からオリゴDNAセンス鎖とアンチセンス鎖を設計して返す
		/// </summary>
		/// <param name="name">標的配列名称</param>
		/// <param name="seq">標的配列（20 bp）</param>
		/// <returns></returns>
		private static Dictionary<string, string> MakeOligoSeq(string name, string seq)
		{
			Dictionary<string, string> res = new Dictionary<string, string>();

			string name_s, name_as;
			name_s = name + "_s";
			name_as = name + "_as";

			// 標的配列は20塩基でなければならない
			char[] seqArray = seq.ToCharArray();
			if (seqArray.Length != 20)
			{
				throw new ArgumentOutOfRangeException();
			}

			// 標的配列の逆相補鎖の作成（全て大文字）
			char[] seqRevCompArray = new char[seqArray.Length];
			char[] seqArrayCap = new char[seqArray.Length];
			for (var i = 0; i < seqArray.Length; i++)
			{
				if ((seqArray[i] == 'G') || (seqArray[i] == 'g'))
				{
					seqRevCompArray[seqArray.Length - i - 1] = 'C';
					seqArrayCap[i] = 'G';
				}
				else if ((seqArray[i] == 'C') || (seqArray[i] == 'c'))
				{
					seqRevCompArray[seqArray.Length - i - 1] = 'G';
					seqArrayCap[i] = 'C';
				}
				else if ((seqArray[i] == 'A') || (seqArray[i] == 'a'))
				{
					seqRevCompArray[seqArray.Length - i - 1] = 'T';
					seqArrayCap[i] = 'A';
				}
				else if ((seqArray[i] == 'T') || (seqArray[i] == 't'))
				{
					seqRevCompArray[seqArray.Length - i - 1] = 'A';
					seqArrayCap[i] = 'T';
				}
				else
				{
					throw new ArgumentOutOfRangeException();
				}
			}

			// センス鎖の設計（追加されるcaccオーバーハングとgは小文字）
			string seqRes, seqRevCompRes;
			if (seqArrayCap[0] != 'G')
			{
				char[] seqArrayG = new char[21];
				seqArrayG[0] = 'g';
				seqArrayCap.CopyTo(seqArrayG, 1);
				seqRes = new string(seqArrayG);
			}
			else
			{
				seqRes = new string(seqArray);
			}
			res.Add(name_s, "cacc" + seqRes);

			// アンチセンス鎖の設計（追加されるaaacオーバーハングとcは小文字）
			if (seqRevCompArray[19] != 'C')
			{
				char[] seqRevCompArrayG = new char[21];
				seqRevCompArray.CopyTo(seqRevCompArrayG, 0);
				seqRevCompArrayG[20] = 'c';
				seqRevCompRes = new string(seqRevCompArrayG);
			}
			else
			{
				seqRevCompRes = new string(seqRevCompArray);
			}
			res.Add(name_as, "aaac" + seqRevCompRes);

			return res;
		}

		private static Dictionary<string, string> AnnealOligos(string name, string seq_s)
        {
			Dictionary<string, string> res = new Dictionary<string, string>();

			string seq_annealed = seq_s + "gttt";
			res.Add(name, seq_annealed);

			return res;
		}

		/// <summary>
		/// 与えられた32ビット整数をリトルエンディアンのバイト列に変換します
		/// </summary>
		/// <param name="input">32ビット整数</param>
		/// <returns>変換後のバイト列</returns>
		private static byte[] Int32ToLittleEndian(int input)
		{
			return ToLittleEndian(BitConverter.GetBytes(input));
		}

		/// <summary>
		/// バイト列のエンディアンを逆転
		/// </summary>
		/// <param name="input">入力バイト列</param>
		/// <returns>エンディアンを逆転したバイト列</returns>
		private static byte[] ToLittleEndian(byte[] input)
        {
			Array.Reverse(input);
			return input;
        }

		/// <summary>
		/// 指定されたDNA配列を二本鎖DNAとしてSnapGeneファイルを作成する
		/// </summary>
		/// <param name="filepath">作成するSnapGeneファイルのフルパス</param>
		/// <param name="seq">DNA配列</param>
		/// <param name="annealed">アニーリング済みならtrue、そうでないならfalse</param>
		private static void MakeSGFile(string filepath, string seq, bool annealed)
        {
			// ヘッダー部分
			byte header_id = 9;
			byte[] header_length = Int32ToLittleEndian(14);
			string header = "SnapGene";
			byte[] header_content = System.Text.Encoding.ASCII.GetBytes(header);
			byte[] header_spacer = { 0x00, 0x01, 0x00, 0x0E, 0x00, 0x0F };

			// 配列部分
			byte seq_id = 0;
			byte[] seq_length = Int32ToLittleEndian(seq.Length + 1);
			byte ds_id = 0b_00000010;
			byte[] seq_content = System.Text.Encoding.ASCII.GetBytes(seq);

			// 末端設定部分
			byte end_id = 8;
			byte[] end_length = Int32ToLittleEndian(263);
			string end_text;
			if (annealed == false)
			{
				end_text = "<AdditionalSequenceProperties><UpstreamStickiness>0</UpstreamStickiness><DownstreamStickiness>0</DownstreamStickiness><UpstreamModification>Unmodified</UpstreamModification><DownstreamModification>Unmodified</DownstreamModification></AdditionalSequenceProperties>";
			}
			else
            {
				end_text = "<AdditionalSequenceProperties><UpstreamStickiness>4</UpstreamStickiness><DownstreamStickiness>4</DownstreamStickiness><UpstreamModification>Unmodified</UpstreamModification><DownstreamModification>Unmodified</DownstreamModification></AdditionalSequenceProperties>";
			}
			byte[] end_content = System.Text.Encoding.ASCII.GetBytes(end_text);

			// ファイルが既に存在するときのエラー処理は未実装
			var bw = new BinaryWriter(new FileStream(filepath, FileMode.CreateNew));
			try
			{
				bw.Write(header_id);
				bw.Write(header_length);
				bw.Write(header_content);
				bw.Write(header_spacer);
				bw.Write(seq_id);
				bw.Write(seq_length);
				bw.Write(ds_id);
				bw.Write(seq_content);
				bw.Write(end_id);
				bw.Write(end_length);
				bw.Write(end_content);
            }
            catch (Exception ex)
            {
				Console.WriteLine("  " + ex.Message);
				Console.WriteLine("  終了するには何かキーを押してください。");
				Console.ReadKey();
				Environment.Exit(0);
			}
            finally
            {
				bw.Close();
            }
        }
	}
}
