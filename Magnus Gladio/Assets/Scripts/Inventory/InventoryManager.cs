using System.Collections;
using System.Collections.Generic;
using System.Security.Cryptography;
using UnityEngine;
using UnityEngine.UI;

public class InventoryManager : MonoBehaviour
{
  public static InventoryManager instance;

  public InventoryScript myBag;
  public GameObject slotGrid;
  public Slot slotPrefab;
  public Text ItemInformation;

    void Awake()
    {
        if (instance != null)
          Destroy(this);
        instance = this;
        myBag.itemList = new List<Item>();
    }

    private void OnEnable()
    {
        RefreshItem();
        instance.ItemInformation.text = "";
    }
    public static void UpdateItemInfo(string itemDescription)
    {
        instance.ItemInformation.text = itemDescription;
    }
    public static void CreateNewItem(Item item)
    {
        Slot newItem = Instantiate(instance.slotPrefab, instance.slotPrefab.transform.position, Quaternion.identity);
        newItem.gameObject.transform.SetParent(instance.slotGrid.transform);
        newItem.slotItem = item;
        newItem.slotImage.sprite = item.itemImage;
       // newItem.slotAmount.text = item.itemAmount.ToString();
       // newItem.slotNumber = instance.myBag.itemList.Count + 1;
        newItem.GetComponent<Button>().onClick.AddListener(newItem.ItemOnclicked);        
    }

    public static void RefreshItem()
    {   

        
        for (int i = 0; i < instance.myBag.itemList.Count; i++)
        {
            CreateNewItem(instance.myBag.itemList[i]);
        }
    }
}
