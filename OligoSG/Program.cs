using System;
using System.IO;
using Microsoft.VisualBasic.FileIO;
using System.Collections.Generic;
using System.Linq;
using System.Xml.Linq;
using System.Text;
using System.Globalization;
using System.Threading.Tasks;

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
			Console.WriteLine("    「標的配列名称+タブ+標的配列（5' -> 3'）+改行」");
			Console.WriteLine("    リストの上下や行間に余計な空白がないように");
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
			for (int i = 0; i < filePath.Length; i++)
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

			// ファイルのディレクトリを取得
			string dirName = Path.GetDirectoryName(filePath[startIndex]);

			// クエリディクショナリ
			var listNameSeq = new Dictionary<string, string>();
			using (TextFieldParser parser = new TextFieldParser(filePath[startIndex], Encoding.ASCII))
			{
				try
				{
					parser.TextFieldType = FieldType.Delimited;
					parser.SetDelimiters("\t");

					while (parser.EndOfData == false)
					{
						string[] data = parser.ReadFields();

						// RESULTSヘッダーがある場合はエラー
						if (data[0] == "RESULTS")
                        {
							Console.WriteLine("  {0}は処理済みのファイルのようです。", filePath[startIndex]);
							Console.WriteLine("  終了するには何かキーを押してください。");
							Console.ReadKey();
							Environment.Exit(0);
						}
						for (int i = 0; i < (data.Length / 2); i++)
						{
							string key, value;
							key = data[i * 2];
							value = data[i * 2 + 1];
							listNameSeq.Add(key, value);
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

			// 結果ディクショナリ
			Dictionary<string, string> oligos = new Dictionary<string, string>();
			Dictionary<string, string> annealed_oligos = new Dictionary<string, string>();
			try
			{
				foreach (KeyValuePair<string, string> pair in listNameSeq)
				{
					Console.WriteLine("  処理中 {0}: {1}", pair.Key, pair.Value);
					try
					{
						Dictionary<string, string> oligo = MakeOligoSeq(pair.Key, pair.Value);
						oligos = oligos.Concat(oligo).ToDictionary(x => x.Key, x => x.Value);
						annealed_oligos = annealed_oligos.Concat(AnnealOligos(pair.Key, oligo[pair.Key + "_s"])).ToDictionary(x => x.Key, x => x.Value);
					}
					catch (ArgumentOutOfRangeException)
					{
						Console.WriteLine("  塩基数は20である必要があります");
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

	class SGFile
    {

    }
}
