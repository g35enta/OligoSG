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

			// クエリディクショナリをインスタンス化
			var listNameSeq = new Dictionary<string, string>();
			// 結果ディクショナリをインスタンス化
			var oligos = new Dictionary<string, string>();
			var annealed_oligos = new Dictionary<string, string>();

			// 入力ファイルを読み込む
			using (TextFieldParser parser = new TextFieldParser(filePath[startIndex], Encoding.ASCII))
			{
				try
				{
					parser.TextFieldType = FieldType.Delimited;
					parser.SetDelimiters("\t");

					while (!parser.EndOfData)
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
						// BbsI単独ベクターの場合
						else if (vector[0] == "pX335" || vector[0] == "pBabe Puro")
						{
							// 標的配列名称と標的配列を取得
							string[] data = parser.ReadFields();
							listNameSeq.Add(data[0], data[1]);

							// コンストラクション用オリゴDNA配列を作成
							Dictionary<string, string> oligo = MakeOligoSeq(data[0], data[1]);

							// 結果ディクショナリを更新
							oligos = oligos.Concat(oligo).ToDictionary(x => x.Key, x => x.Value);
							annealed_oligos = annealed_oligos.Concat(AnnealOligos(data[0], oligo[data[0] + "_s"])).ToDictionary(x => x.Key, x => x.Value);

							// ベクターファイルを読み込むリスト
							List<SGSegment> vec = new List<SGSegment>();

							// 配列挿入部位を読み込むリスト
							List<byte> pattern = new List<byte>();

							// 5'末端にGが付加されているか否かで分岐
							if (oligo[data[0] + "_s"].Length == 24)
							{
								vec = ReadSegments(exeDirName + "\\" + vector[0] + "_N20.dna");
								pattern.AddRange(System.Text.Encoding.ASCII.GetBytes("NNNNNNNNNNNNNNNNNNNN").ToList<byte>()); // N20
							}
							else if (oligo[data[0] + "_s"].Length == 25)
							{
								vec = ReadSegments(exeDirName + "\\" + vector[0] + "_N21.dna");
								pattern.AddRange(System.Text.Encoding.ASCII.GetBytes("NNNNNNNNNNNNNNNNNNNNN").ToList<byte>()); // N21
							}
							else { }

							// 配列挿入部位を標的配列で置換してSnapGeneセグメントをリストに格納
							List<byte> seq_original = new List<byte>();
							foreach (SGSegment seg in vec)
							{
								if (seg.SGIdentifier == 0)
								{
									seq_original = seg.SGContent.ToList<byte>();

									// センス鎖から5'オーバーハングを除去した配列を取得
									string seq_s_nooverhangs = oligo[data[0] + "_s"].Remove(0, 4);
									byte[] seqBytes = System.Text.Encoding.ASCII.GetBytes(seq_s_nooverhangs);

									// ベクター配列内のNNN…NNNをオリゴ配列で置換したものを取得
									List<byte> seq_replaced = ReplaceList(seq_original, pattern, seqBytes.ToList<byte>());

									seg.SGContent = seq_replaced.ToArray();
								}
							}

							// SnapGeneファイルとして保存
							string path = dirName + "\\" + data[0] + "_" + vector[0] + ".dna";
							SaveSGFiles(vec, path);
							Console.WriteLine("  ベクターマップを作成しました：{0}", path);
						}
						// BbsI/BsmBI両方あるベクターの場合
						else if (vector[0] == "pX459dual D10A" || vector[0] == "pX459dual EGFP D10A")
                        {
							// 標的配列名称と標的配列を二組取得
							string[] data1 = parser.ReadFields();
							listNameSeq.Add(data1[0], data1[1]);
							string[] data2 = parser.ReadFields();
							listNameSeq.Add(data2[0], data2[1]);

							// コンストラクション用オリゴDNA配列を作成
							Dictionary<string, string> oligo1 = MakeOligoSeq(data1[0], data1[1]);
							Dictionary<string, string> oligo2 = MakeOligoSeq(data2[0], data2[1]);

							// 結果ディクショナリを更新
							oligos = oligos.Concat(oligo1).ToDictionary(x => x.Key, x => x.Value);
							oligos = oligos.Concat(oligo2).ToDictionary(x => x.Key, x => x.Value);
							annealed_oligos = annealed_oligos.Concat(AnnealOligos(data1[0], oligo1[data1[0] + "_s"])).ToDictionary(x => x.Key, x => x.Value);
							annealed_oligos = annealed_oligos.Concat(AnnealOligos(data2[0], oligo2[data2[0] + "_s"])).ToDictionary(x => x.Key, x => x.Value);

							// ベクターファイルを読み込むリスト
							List<SGSegment> vec = new List<SGSegment>();

							// 配列挿入部位を読み込むリスト
							List<byte> pattern_BbsI = new List<byte>();
							List<byte> pattern_BsmBI = new List<byte>();

							// 5'末端にGが付加されているか否かで分岐
							if (oligo1[data1[0] + "_s"].Length == 24) // BbsI
							{
								if (oligo2[data2[0] + "_s"].Length == 24) // BsmBI
								{
									vec = ReadSegments(exeDirName + "\\" + vector[0] + "_BsmBIM20_BbsIN20.dna");
									pattern_BbsI.AddRange(System.Text.Encoding.ASCII.GetBytes("NNNNNNNNNNNNNNNNNNNN").ToList<byte>()); // N20
									pattern_BsmBI.AddRange(System.Text.Encoding.ASCII.GetBytes("MMMMMMMMMMMMMMMMMMMM").ToList<byte>()); // M20
								}
								else if (oligo2[data2[0] + "_s"].Length == 25) // BsmBI
                                {
									vec = ReadSegments(exeDirName + "\\" + vector[0] + "_BsmBIM21_BbsIN20.dna");
									pattern_BbsI.AddRange(System.Text.Encoding.ASCII.GetBytes("NNNNNNNNNNNNNNNNNNNN").ToList<byte>()); // N20
									pattern_BsmBI.AddRange(System.Text.Encoding.ASCII.GetBytes("MMMMMMMMMMMMMMMMMMMMM").ToList<byte>()); // M21
								}
								else { }
							}
							else if (oligo1[data1[0] + "_s"].Length == 25) // BbsI
							{
								if (oligo2[data2[0] + "_s"].Length == 24) // BsmBI
								{
									vec = ReadSegments(exeDirName + "\\" + vector[0] + "_BsmBIM20_BbsIN21.dna");
									pattern_BbsI.AddRange(System.Text.Encoding.ASCII.GetBytes("NNNNNNNNNNNNNNNNNNNNN").ToList<byte>()); // N21
									pattern_BsmBI.AddRange(System.Text.Encoding.ASCII.GetBytes("MMMMMMMMMMMMMMMMMMMM").ToList<byte>()); // M20
								}
								else if (oligo2[data2[0] + "_s"].Length == 25) // BsmBI
								{
									vec = ReadSegments(exeDirName + "\\" + vector[0] + "_BsmBIM21_BbsIN21.dna");
									pattern_BbsI.AddRange(System.Text.Encoding.ASCII.GetBytes("NNNNNNNNNNNNNNNNNNNNN").ToList<byte>()); // N21
									pattern_BsmBI.AddRange(System.Text.Encoding.ASCII.GetBytes("MMMMMMMMMMMMMMMMMMMMM").ToList<byte>()); // M21
								}
								else { }
							}
							else { }

							// 配列挿入部位を標的配列で置換してSnapGeneセグメントをリストに格納
							List<byte> seq_original = new List<byte>();
							foreach (SGSegment seg in vec)
							{
								if (seg.SGIdentifier == 0)
								{
									seq_original = seg.SGContent.ToList<byte>();

									// センス鎖から5'オーバーハングを除去した配列を取得
									string seq_s_nooverhangs1 = oligo1[data1[0] + "_s"].Remove(0, 4);
									string seq_s_nooverhangs2 = oligo2[data2[0] + "_s"].Remove(0, 4);
									byte[] seqBytes_BbsI = System.Text.Encoding.ASCII.GetBytes(seq_s_nooverhangs1);
									byte[] seqBytes_BsmBI = System.Text.Encoding.ASCII.GetBytes(seq_s_nooverhangs2);

									// ベクター配列内の挿入部位をオリゴ配列で置換したものを取得
									List<byte> seq_replaced_BbsI = ReplaceList(seq_original, pattern_BbsI, seqBytes_BbsI.ToList<byte>());
									List<byte> seq_replaced_BbsI_BsmBI = ReplaceList(seq_replaced_BbsI, pattern_BsmBI, seqBytes_BsmBI.ToList<byte>());

									seg.SGContent = seq_replaced_BbsI_BsmBI.ToArray();
								}
							}

							// SnapGeneファイルとして保存
							string path = dirName + "\\" + data1[0] + "_" + data2[0] + "_" + vector[0] + ".dna";
							SaveSGFiles(vec, path);
							Console.WriteLine("  ベクターマップを作成しました：{0}", path);
						}
						else
                        {
							Console.WriteLine("  入力ファイルの形式が不正です。");
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

			byte[] HEADER_BYTE = { 0x53, 0x6E, 0x61, 0x70, 0x47, 0x65, 0x6E, 0x65, 0x00, 0x01, 0x00, 0x0E, 0x00, 0x10 };

			// それぞれのオリゴDNA配列SnapGeneファイル作成
			foreach (KeyValuePair<string, string> kv in oligos)
			{
				try
				{
					string path = dirName + "\\" + kv.Key + ".dna";
					byte[] seq_byte = System.Text.Encoding.ASCII.GetBytes("2" + kv.Value);
                    byte[] no_stickiness_byte = System.Text.Encoding.ASCII.GetBytes("<AdditionalSequenceProperties><UpstreamStickiness>0</UpstreamStickiness><DownstreamStickiness>0</DownstreamStickiness><UpstreamModification>Unmodified</UpstreamModification><DownstreamModification>Unmodified</DownstreamModification></AdditionalSequenceProperties>");
					SGSegment header = new SGSegment(0x09, 14, HEADER_BYTE);
					SGSegment seq = new SGSegment(0x00, kv.Value.Length + 1, seq_byte);
					SGSegment end = new SGSegment(0x08, 263, no_stickiness_byte);
                    List<SGSegment> oligoSGSegs = new List<SGSegment>
                    {
                        header,
                        seq,
                        end
                    };
                    SaveSGFiles(oligoSGSegs, path);

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
					byte[] seq_byte = System.Text.Encoding.ASCII.GetBytes("2" + kv.Value);
					byte[] no_stickiness_byte = System.Text.Encoding.ASCII.GetBytes("<AdditionalSequenceProperties><UpstreamStickiness>4</UpstreamStickiness><DownstreamStickiness>4</DownstreamStickiness><UpstreamModification>Unmodified</UpstreamModification><DownstreamModification>Unmodified</DownstreamModification></AdditionalSequenceProperties>");
					SGSegment header = new SGSegment(0x09, 14, HEADER_BYTE);
					SGSegment seq = new SGSegment(0x00, kv.Value.Length + 1, seq_byte);
					SGSegment end = new SGSegment(0x08, 263, no_stickiness_byte);
                    List<SGSegment> annealedoligoSGSegs = new List<SGSegment>
                    {
                        header,
                        seq,
                        end
                    };
                    SaveSGFiles(annealedoligoSGSegs, path);

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

			Console.WriteLine("  全ての処理が完了しました。");
			Console.WriteLine("  終了するには何かキーを押してください。");
			Console.ReadKey();
			Environment.Exit(0);
		}
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
		/// DNA配列フラグを文字列に変換して返します
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

		/// <summary>
		/// 与えられたバイト型配列のなかで与えられたパターンに一致する部分を与えられた配列で置換します
		/// </summary>
		/// <param name="input">置換対象の配列全体</param>
		/// <param name="pattern">置換される配列</param>
		/// <param name="replacement">置き換わる配列</param>
		/// <returns>置換後の配列全体</returns>
		public static byte[] Replace(byte[] input, byte[] pattern, byte[] replacement)
		{
			if (pattern.Length == 0)
			{
				return input;
			}

			List<byte> result = new List<byte>();

			// 検索開始地点
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

				// pattern全長と一致する部分があった場合のみfoundMatchがtrueに維持される
				// 一致する部分があった場合には、replacementを結果に追加して、
				// 検索開始地点をpatternのぶん飛ばす
				if (foundMatch)
				{
					result.AddRange(replacement);
					i += pattern.Length - 1; // 次のループで加算されるぶん1少ない
				}
				// 検索開始地点からpatternの文字数だけ比較してすべて一致しなかった場合には、
				// 検索開始地点の要素をそのまま結果に追加する
				else
				{
					result.Add(input[i]);
				}
			}

			// inpurの終わりのpatternの長さ未満の配列はそのままresultに追加する
			for (/**/; i < input.Length; i++)
			{
				result.Add(input[i]);
			}

			return result.ToArray();
		}

		/// <summary>
		/// 与えられたバイト型配列のなかで与えられたパターンに一致する部分を与えられた配列で置換します
		/// </summary>
		/// <param name="input">置換対象の配列全体</param>
		/// <param name="pattern">置換される配列</param>
		/// <param name="replacement">置き換わる配列</param>
		/// <returns>置換後の配列全体</returns>
		public static List<byte> ReplaceList(List<byte> input, List<byte> pattern, List<byte> replacement)
		{
			if (pattern.Count == 0)
			{
				return input;
			}

			List<byte> result = new List<byte>();

			// 検索開始地点
			int i;
			for (i = 0; i <= input.Count - pattern.Count; i++)
			{
				bool foundMatch = true;
				for (int j = 0; j < pattern.Count; j++)
				{
					if (input[i + j] != pattern[j])
					{
						foundMatch = false;
						break;
					}
				}

				// pattern全長と一致する部分があった場合のみfoundMatchがtrueに維持される
				// 一致する部分があった場合には、replacementを結果に追加して、
				// 検索開始地点をpatternのぶん飛ばす
				if (foundMatch)
				{
					result.AddRange(replacement);
					i += pattern.Count - 1; // 次のループで加算されるぶん1少ない
				}
				// 検索開始地点からpatternの文字数だけ比較してすべて一致しなかった場合には、
				// 検索開始地点の要素をそのまま結果に追加する
				else
				{
					result.Add(input[i]);
				}
			}

			// inpurの終わりのpatternの長さ未満の配列はそのままresultに追加する
			for (/**/; i < input.Count; i++)
			{
				result.Add(input[i]);
			}

			return result;
		}
		/// <summary>
		/// 与えられた配列からオリゴDNAセンス鎖とアンチセンス鎖を設計して返す
		/// </summary>
		/// <param name="name">標的配列名称</param>
		/// <param name="seq">標的配列（20 bp）</param>
		/// <returns></returns>
		public static Dictionary<string, string> MakeOligoSeq(string name, string seq)
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

		/// <summary>
		/// 与えられたセンス鎖配列からアニーリング状態の配列を返す
		/// </summary>
		/// <param name="name">標的配列名称</param>
		/// <param name="seq_s">アニーリング状態の配列</param>
		/// <returns></returns>
		public static Dictionary<string, string> AnnealOligos(string name, string seq_s)
        {
			Dictionary<string, string> res = new Dictionary<string, string>();

			string seq_annealed = seq_s + "gttt";
			res.Add(name, seq_annealed);

			return res;
		}
	}
}
